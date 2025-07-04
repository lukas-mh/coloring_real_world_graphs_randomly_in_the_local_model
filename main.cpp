#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>
#include <set>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <atomic>
#include <map>
#include <cassert>
#include <unistd.h>
#include <variant>

#include "girgs/source/girgs/include/girgs/Generator.h"

#define NUMBER_OF_VERTICES 100
#define NUMBER_OF_COLORS 50
#define AVERAGE_DEGREE 10
#define PI 3.1415926
#define RANDOM_GRAPH_SCALE 10000

bool debug = false;
bool extra_detail_output = debug && true;

using namespace std;

class Graph;

vector<variant<ulong, double, std::string, list<int>>> color_randomly(double epsilon_, std::string factoring_, int random_seed_, Graph graph);

class Graph
{
    vector<pair<int, int>> edges;
    int n;
    int avg_degree;
    double alpha;
    int graph_seed;

public:
    Graph(vector<pair<int, int>> edges_, int n_, int avg_degree_, double alpha_, int graph_seed_) : edges(edges_), n(n_), avg_degree(avg_degree_), alpha(alpha_), graph_seed(graph_seed_) {};

    vector<variant<ulong, double, std::string, list<int>>> color(double epsilon, std::string factoring, int random_seed)
    {
        return color_randomly(epsilon, factoring, random_seed, *this);
    }

    vector<pair<int, int>> const get_edges() { return edges; }
    int const get_n() { return n; }
    int const get_average_degree() { return avg_degree; }
    double const get_alpha() { return alpha; }
    int const get_graph_seed() { return graph_seed; }
};

enum ColorState
{
    UNCOLORED,
    COLORED,
    CANDIDATE
};

struct V
{
    ulong index;
    ulong color;
    ColorState color_state = {UNCOLORED};
    unordered_set<V *> neighbours;
    unordered_set<ulong> neighbours_colors;
    vector<ulong> neighbours_colors_v;

    V(ulong _index, ulong delta, ulong _color)
    {
        index = _index;
        color = _color;
        neighbours = unordered_set<V *>(delta);
        neighbours_colors = unordered_set<ulong>(delta);
        neighbours_colors_v = vector<ulong>(delta);
    }
    V(ulong _index, ulong delta)
    {
        index = _index;
        neighbours = unordered_set<V *>(delta);
        neighbours_colors = unordered_set<ulong>(delta);
        neighbours_colors_v = vector<ulong>(delta);
    }
    V(ulong _index) { index = _index; }
    V() { index = 0; }
};

void generate_conflict_graph(list<V> &vertices, list<V *> &vertices_p, vector<ulong> &deltas, std::vector<std::pair<int, int>> &edges, std::vector<std::pair<int, int>> &conflicted_edges, vector<ulong> &v_colors, int n)
{
    auto t_conflicted_indices = chrono::high_resolution_clock::now();
    set<int> conflicted_indices;
    vector<bool> index_is_conflicted(n, false);
    for (auto &conflicted_edge : conflicted_edges)
    {
        index_is_conflicted[conflicted_edge.first] = true;
        index_is_conflicted[conflicted_edge.second] = true;
        conflicted_indices.emplace(conflicted_edge.first);
        conflicted_indices.emplace(conflicted_edge.second);
    }
    auto t_conflicted_indices_end = chrono::high_resolution_clock::now();

    auto t_conflict_map = chrono::high_resolution_clock::now();
    vertices.clear();
    map<int, V *> conflicted_vertices_map;

    for (int idx : conflicted_indices)
    {
        vertices.emplace_back(V(idx, deltas.at(idx), v_colors.at(idx)));
        conflicted_vertices_map.emplace(idx, &vertices.back());
        vertices_p.emplace_back(&vertices.back());
    }
    auto t_conflict_map_end = chrono::high_resolution_clock::now();

    map<int, int> conflict_deltas;
    if (extra_detail_output)
    {
        for (auto ci : conflicted_indices)
        {
            conflict_deltas.emplace(ci, 0);
        }
    }

    auto t_fill_edges = chrono::high_resolution_clock::now();
    for (auto const &edge : edges) // ALL edges
    {
        // auto first = conflicted_vertices_map.find(edge.first); // first vertex of edge
        // auto second = conflicted_vertices_map.find(edge.second); // second vertex of edge
        if (index_is_conflicted[edge.first]) // in case it is a conflicted vertex (could be found in exactly this map)
        {
            if (index_is_conflicted[edge.second])
            { // if second one is confliced as well it's considered a neighbour
                conflicted_vertices_map.find(edge.first)->second->neighbours.emplace(conflicted_vertices_map.find(edge.second)->second);
            }
            else
            { // otherwise just a forbidden color
                conflicted_vertices_map.find(edge.first)->second->neighbours_colors_v.push_back(v_colors.at(edge.second));
            }
            if (extra_detail_output)
            {
                conflict_deltas.at(edge.first) = conflict_deltas.at(edge.first) + 1;
            }
        }
        if (index_is_conflicted[edge.second]) // see first one
        {
            if (index_is_conflicted[edge.first])
            {
                conflicted_vertices_map.find(edge.second)->second->neighbours.emplace(conflicted_vertices_map.find(edge.first)->second);
            }
            else
            {
                conflicted_vertices_map.find(edge.second)->second->neighbours_colors_v.push_back(v_colors.at(edge.first));
            }
            if (extra_detail_output)
            {
                conflict_deltas.at(edge.second) = conflict_deltas.at(edge.second) + 1;
            }
        }
    }

    for (auto &v : vertices)
    {
        unordered_set<ulong> neighbours_colors_set(v.neighbours_colors_v.begin(), v.neighbours_colors_v.end());
        v.neighbours_colors_v.clear();
        v.neighbours_colors = neighbours_colors_set;
    }
    auto t_fill_edges_end = chrono::high_resolution_clock::now();

    return;
}

Graph generate_graph(int n, int graph_seed, int d = 2, double ple = 2.5, double alpha = std::numeric_limits<double>::infinity(), double deg = 10, int wseed = 12, int pseed = 130, int sseed = 1400, int threads = 10)
{
    auto weights = girgs::generateWeights(n, ple, wseed, true /* threads > 1 */);
    auto positions = girgs::generatePositions(n, d, pseed);
    auto scaling = girgs::scaleWeights(weights, deg, d, alpha);
    vector<pair<int, int>> edges = girgs::generateEdges(weights, positions, alpha, sseed);

    return Graph(edges, n, deg, alpha, graph_seed);
}

vector<ulong> get_delta_edges(vector<std::pair<int, int>> &edges, int n)
{
    vector<ulong> delta_counter(n);
    for (std::pair<int, int> &edge : edges)
    {
        delta_counter.at(edge.first)++;
        delta_counter.at(edge.second)++;
    }
    return delta_counter;
}

vector<variant<ulong, double, std::string, list<int>>> color_randomly(double epsilon_, std::string factoring_, int random_seed_, Graph graph)
{
    auto t_total = chrono::high_resolution_clock::now();

    // parsing arguments
    int number_of_vertices = graph.get_n();
    double epsilon = epsilon_;
    int graph_seed = graph.get_graph_seed();
    int random_seed = random_seed_;
    int avg_degree = graph.get_average_degree();
    double alpha = graph.get_alpha();
    std::string factoring = factoring_;

    // graph variables
    list<V> vertices;
    list<V *> vertices_p;
    std::vector<std::pair<int, int>> edges = graph.get_edges();
    vector<ulong> v_colors(number_of_vertices);

    // calculating delta and c (number of colors)
    vector<ulong> deltas = get_delta_edges(edges, number_of_vertices);
    ulong delta = *max_element(deltas.begin(), deltas.end());
    ulong c;

    if (factoring == "mul")
    {
        c = (int)((double)delta * epsilon);
    }
    else if (factoring == "exp")
    {
        c = (int)pow(delta, epsilon);
    }
    else
    {
        exit(1);
        list l = {0};
        vector<variant<ulong, double, std::string, list<int>>> ret = {(ulong)number_of_vertices, (ulong)delta, (double)epsilon, factoring, (ulong)0, (ulong)0, (ulong)0, l, (ulong)avg_degree, (double)alpha, (ulong)graph_seed, (ulong)random_seed, (ulong)c};
        return ret;
    }

    if (c < 1)
    {
        list l = {0};
        vector<variant<ulong, double, std::string, list<int>>> ret = {(ulong)number_of_vertices, (ulong)delta, (double)epsilon, factoring, (ulong)0, (ulong)0, (ulong)0, l, (ulong)avg_degree, (double)alpha, (ulong)graph_seed, (ulong)random_seed, (ulong)c};
        return ret;
    }

    int random_nr = random_seed;
    srand(random_nr);

    list<int> conflicts_progression;

    // First round
    for (ulong &v : v_colors)
    {
        v = rand() % c; // very minor caveat: this leads to e.g. (c=10_000) the first 3647 numbers beeing 0.000465662% more likely than later ones
    }
    vector<pair<int, int>> conflicted_edges;
    for (auto const &edge : edges)
    {
        if (v_colors.at(edge.first) == v_colors.at(edge.second))
        {
            conflicted_edges.emplace_back(edge);
        }
    }
    bool coloring_complete = conflicted_edges.empty();
    int runtime_counter = 1;
    if (coloring_complete)
    {
        int is_success = coloring_complete;
        conflicts_progression.emplace_back(0);
        vector<variant<ulong, double, std::string, list<int>>> ret = {(ulong)number_of_vertices, (ulong)delta, (double)epsilon, factoring, (ulong)runtime_counter, (ulong)is_success, (ulong)(conflicts_progression.back()), conflicts_progression, (ulong)avg_degree, (double)alpha, (ulong)graph_seed, (ulong)random_nr, (ulong)c};
        return ret;
    }

    generate_conflict_graph(vertices, vertices_p, deltas, edges, conflicted_edges, v_colors, number_of_vertices);

    conflicts_progression.emplace_back(vertices.size());

    // Big While Loop
    while (!coloring_complete && runtime_counter < 1000)
    {
        int max_roll_counter = 0;
        for (auto vertex : vertices_p)
        {
            int roll_counter = 0;
            int rand_color;
            do
            {
                rand_color = rand() % c; // very minor caveat: this leads to e.g. (c=10_000) the first 3647 numbers beeing 0.000465662% more likely than later ones
                roll_counter++;
            } while (vertex->neighbours_colors.count(rand_color) > 0 && roll_counter < 10000000 && vertex->neighbours_colors.size() < c);
            if (roll_counter >= 10000000)
            {
                int is_success = 0;
                vector<variant<ulong, double, std::string, list<int>>> ret = {(ulong)number_of_vertices, (ulong)delta, (double)epsilon, factoring, (ulong)runtime_counter, (ulong)is_success, (ulong)(conflicts_progression.back()), conflicts_progression, (ulong)avg_degree, (double)alpha, (ulong)graph_seed, (ulong)random_nr, (ulong)c};
                return ret;
            }
            vertex->color = rand_color;
            vertex->color_state = CANDIDATE;
            max_roll_counter = roll_counter > max_roll_counter ? roll_counter : max_roll_counter;
        }

        list<V *> vertices_p_new = {};
        vertices_p_new.clear();
        atomic<int> conflicts_counter(0);
        for (auto vertex_it = vertices_p.begin(); vertex_it != vertices_p.end();)
        {
            if ((*vertex_it)->color_state == CANDIDATE)
            {
                for (V *neighbour : (*vertex_it)->neighbours)
                {
                    if ((*vertex_it)->color == neighbour->color)
                    {
                        neighbour->color_state = UNCOLORED;
                        (*vertex_it)->color_state = UNCOLORED;
                        break;
                    }
                }
            }
            if ((*vertex_it)->color_state == UNCOLORED)
            {
                vertices_p_new.emplace_back((*vertex_it));
                conflicts_counter += 1;
            }
            else
            {
                (*vertex_it)->color_state = COLORED;
            }
            vertex_it++;
        }
        if (vertices_p_new.size() > 0)
        {
            coloring_complete = false;
        }
        else
        {
            coloring_complete = true;
        }
        conflicts_progression.emplace_back(conflicts_counter);

        for (auto vertex : vertices_p_new)
        {
            auto vertex_neighbours = vertex->neighbours;
            for (V *neighbour : vertex_neighbours)
            {
                if (neighbour->color_state == COLORED)
                {
                    vertex->neighbours_colors.insert(neighbour->color);
                    vertex->neighbours.erase(neighbour);
                }
            }
        }
        vertices_p = vertices_p_new;
        auto t_for_3_end = chrono::high_resolution_clock::now();
        runtime_counter++;
    }

    int is_success = coloring_complete;
    vector<variant<ulong, double, std::string, list<int>>> ret = {(ulong)number_of_vertices, (ulong)delta, (double)epsilon, factoring, (ulong)runtime_counter, (ulong)is_success, (ulong)(conflicts_progression.back()), conflicts_progression, (ulong)avg_degree, (double)alpha, (ulong)graph_seed, (ulong)random_nr, (ulong)c};
    return ret;
}
