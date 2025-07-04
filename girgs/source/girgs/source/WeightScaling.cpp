#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>
#include <cmath>
#include <vector>
#include <limits>

#include <girgs/WeightScaling.h>

namespace girgs {

// helper for scale weights
static double exponentialSearch(const std::function<double(double)> &f, double desiredValue, double accuracy = 0.02, double lower = 1.0,
                                double upper = 2.0) {

    // scale interval up if necessary
    while (f(upper) < desiredValue) {
        lower = upper;
        upper *= 2;
    }

    // scale interval down if necessary
    while (f(lower) > desiredValue) {
        upper = lower;
        lower /= 2;
    }

    // do binary search
    auto mid = f((upper + lower) / 2);
    while (std::abs(mid - desiredValue) > accuracy) {
        if (mid < desiredValue)
            lower = (upper + lower) / 2;
        else
            upper = (upper + lower) / 2;
        mid = f((upper + lower) / 2);
    }

    return (upper + lower) / 2;
}


/// Safe some time by not sorting the complete weights vector, but only the largest elements
class LazySorter {
public:
    using iterator = typename std::vector<double>::iterator;

    LazySorter(std::vector<double>& sweights):
        m_sweights(sweights), m_sorted_end(sweights.begin())
    {}

    /// Ensures that
    iterator sort_downto(double thresh) {
        if (m_lower > thresh) {
            // we need more points, make sure we have more!
            assert(m_sorted_end != m_sweights.end());

            // move all not yet sorted elements larger than the threshold directly next to the sorted ones
            auto new_ps_end = std::partition(m_sorted_end, m_sweights.end(), [=](double x) { return x >= thresh; });
            if (new_ps_end == m_sorted_end) return m_sorted_end;

            std::sort(m_sorted_end, new_ps_end, std::greater<double>());
            assert(std::is_sorted(m_sweights.begin(), m_sorted_end, std::greater<double>()));

            m_sorted_end = new_ps_end;
            m_lower = (m_sorted_end == m_sweights.end()) ? std::numeric_limits<double>::min() : thresh;
            return m_sorted_end;

        } else {
            // the splitter is within our already sorted segment
            return std::lower_bound(m_sweights.begin(), m_sorted_end, thresh,
                                    [&] (double x, double thresh) {return x > thresh;});

        }
    };


private:
    std::vector<double>& m_sweights;
    iterator m_sorted_end;
    double m_lower{std::numeric_limits<double>::max()};
};

// helper for scale weights
double estimateWeightScalingThreshold(const std::vector<double> &weights, double desiredAvgDegree, int dimension) {
    std::vector<double> sweights(weights.size());
    LazySorter lazy_sorter(sweights);

    // compute some constant stuff
    const auto n = weights.size();
    auto max_weight = 0.0;
    auto W = 0.0, sq_W = 0.0;
    {
#ifndef _MSC_VER
        #pragma omp parallel for reduction(+:W, sq_W), reduction(max: max_weight)
#endif
        for (int i = 0; i < n; ++i) {
            const auto each = weights[i];
            sweights[i] = each; // copy to sweights

            W += each;
            sq_W += each * each;
            max_weight = std::max(max_weight, each);
        }
    }

    // my function to do the exponential search on
    auto f = [=, &sweights, &lazy_sorter](double c) {
        
        // compute overestimation
        const auto pow2c = std::pow(2.0 * c, dimension);
		const auto overestimation = pow2c * (W - sq_W / W);

		// compute rich club
		const auto richclub_end = lazy_sorter.sort_downto(W / pow2c / max_weight);
		assert(richclub_end <= sweights.end());

        // subtract error
        auto error = 0.0;

        auto w2sum = 0.0;
        auto w2 = sweights.cbegin();
        for(auto w1 = richclub_end; w1 > sweights.begin(); ) {
			--w1; // msvc throws when decrementing beyond begin iterator
            const auto fac = *w1 / W;
            const auto my_thres = 1.0 / fac / pow2c;

            for(; w2 < richclub_end && *w2 >= my_thres; w2sum += *(w2++));

            /**
              * sum_{k < j, k != i}{ std::pow(2*c,dimension)*(w1*w_k/W)-1.0 }
              * = sum_{k < j, k != i}{ std::pow(2*c,dimension)*(w1*w_k/W) } - j
              * = sum_{k < j}{ std::pow(2*c,dimension)*(w1*w_k/W) } - j - (0 if j < i else std::pow(2*c,dimension)*(w1*w_i/W) - 1)
              * = pow2c * w1/W * (sum_j{w_j}) - j - (0 if j < i else pow2c*(w1/W)*w_i - 1)
              */

            error += w2sum * pow2c * fac - (std::distance(sweights.cbegin(), w2));

            if (w2 >= w1) {
                // we have to subtract the self-contribution of x == y
                error -= pow2c * fac * (*w1) - 1.0;
            }
        }

        return (overestimation - error) / n;
    };

    // do exponential search on expected average degree function
    auto estimated_c = exponentialSearch(f, desiredAvgDegree);

    /*
     * edge iff dist < c(wi*wj/W)^(1/d)
     *
     * c(wi*wj/W)^(1/d)
     * = ( c^d * wi*wj/W )^(1/d)
     * = ( c^d*wi * c^d*wi / (c^d*W) )^(1/d)
     *
     * so we can just scale all weights by c^d
     */
    return pow(estimated_c, dimension); // return scaling
}


// helper for scale weights
double estimateWeightScaling(const std::vector<double> &weights, double desiredAvgDegree, int dimension, double alpha) {
    assert(alpha != 1.0); // somehow breaks for alpha 1.0

    // compute some constant stuff
    auto W = std::accumulate(weights.begin(), weights.end(), 0.0);
    auto sum_sq_w = 0.0; // sum_{v\in V} (w_v^2/W)
    auto sum_w_a = 0.0; // sum_{v\in V} (w_v  /W)^\alpha
    auto sum_sq_w_a = 0.0; // sum_{v\in V} (w_v^2/W)^\alpha

    //   sum_{u\in V} sum_{v\in V} (wu*wv/W)^\alpha
    // = sum_{u\in V} sum_{v\in V} wu^\alpha * (wv/W)^\alpha
    // = sum_{u\in V} wu^\alpha sum_{v\in V} (wv/W)^\alpha
    auto sum_wwW_a = 0.0;
    auto max_w = 0.;
    const auto n = static_cast<int>(weights.size());

    // this loop causes >= 70% of runtime
    std::vector<double> sweights(n);
    LazySorter lazy_sorter(sweights);
    {
#ifndef _MSC_VER
        #pragma omp parallel for reduction(+:sum_sq_w, sum_w_a, sum_sq_w_a, sum_wwW_a), reduction(max:max_w)
#endif
        for (int i = 0; i < n; ++i) {
            const auto each = weights[i];
            sweights[i] = each; // copy in parallel

            const auto each_W = each / W;
            const auto pow_each = pow(each, alpha);
            const auto pow_each_W = pow(each / W, alpha);

            sum_sq_w += each * each_W;
            sum_wwW_a += pow_each;
            sum_w_a += pow_each_W;
            sum_sq_w_a += pow_each * pow_each_W;
            max_w = std::max(each, max_w);
        }
    }
    sum_wwW_a *= sum_w_a;
    const auto max_w_W = max_w / W;

    const auto factor1 = (W - sum_sq_w) * (1 + 1 / (alpha - 1)) * (1 << dimension);
    const auto factor2 = pow(2, alpha * dimension) / (alpha - 1) * (sum_wwW_a - sum_sq_w_a);

    const auto W_alpha = std::pow(W, alpha);

    // pow_weights[i] = std::pow(sweights[i], alpha) for all elements of the rich club
    std::vector<double> pow_weights;

    auto f = [&] (double c) {
        const auto long_and_short_with_error = pow(c, 1 / alpha) * factor1 - c * factor2;

        const auto rich_thresh = std::exp(dimension * std::log(0.5 / std::pow(c, 1.0 / alpha / dimension)) - log(max_w_W));
        const auto richclub_end = lazy_sorter.sort_downto(rich_thresh);
        const auto num_richclub = static_cast<int>(std::distance(sweights.begin(), richclub_end));

        if (!num_richclub)
            return long_and_short_with_error / n;

        assert(num_richclub <= n);

        // precompute new pows in case the richclub grew
        pow_weights.reserve(n / num_richclub > 2 ? num_richclub : n);
        while(pow_weights.size() < num_richclub) {
            pow_weights.push_back(std::pow(sweights[pow_weights.size()], alpha));
        }

        // get error for long and short edges
        const auto thresh = std::exp( (std::log(0.5) * dimension - std::log(c) / alpha) );

        int i2 = 0;
        auto w2_sum       = 0.0;
        long double w2_alpha_sum = 0.0;

        auto w_terms = 0.0;
        auto w_alpha_terms = 0.0;
        size_t num_terms = 0;

        for(auto i1 = num_richclub; i1-- > 0; ) {
            const auto my_thres = thresh * W / sweights[i1];
            for (; i2 < num_richclub && sweights[i2] >= my_thres; ++i2) {
                w2_sum       += sweights[i2];
                w2_alpha_sum += pow_weights[i2];
            }

            num_terms     += i2;
            w_terms       += sweights[i1]    * w2_sum       / W;
            w_alpha_terms += pow_weights[i1] * w2_alpha_sum / W_alpha;

            if (i2 >= i1) {
                num_terms--;
                w_terms       -= sweights[i1]    * sweights[i1]    / W;
                w_alpha_terms -= pow_weights[i1] * pow_weights[i1] / W_alpha;
            }
        }

        auto short_error = (1 << dimension) * pow(c, 1 / alpha) * w_terms - num_terms;
        auto long_error = c * dimension * (1 << dimension) / (dimension - alpha * dimension) *
            (std::pow(0.5, dimension - alpha * dimension) * w_alpha_terms -
                std::pow(c, 1.0 / alpha - 1.0) * w_terms);

        return (long_and_short_with_error - short_error - long_error) / n;
    };

    // do exponential search on avg_degree function
    auto estimated_c = exponentialSearch(f, desiredAvgDegree);

    /*
     * Pr(edge) = Pr(c * 1/dist^ad * (wi*wj/W)^a )
     *
     * c * (wi*wj/W)^a
     * = (c^{1/a} wi*wj/W)^a
     * = (c^{1/a}wi* (c^{1/a}wj / (c^{1/a}W)^a
     *
     * so we can just scale all weights by (c^{1/a}
     */
    return pow(estimated_c, 1 / alpha); // return scaling
}

} // namespace girgs
