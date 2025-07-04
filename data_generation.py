from build.coloringgirgs import *
import pandas as pd

data_labels = ["n", "maxdelta", "epsilon", "factoring", "runtime", "success", "conflicts_left", "conflict_progression", "avg_degree", "alpha", "graph_seed", "random_seed", "colors"]


def color(args):
    g, e, factoring, random_seed = args
    return g.color(e, factoring, random_seed)


def data_gen(graph_seed_range, random_seed_range, factoring, epsilons, vertices, name=""):
    data = []
    for n in vertices:
        for graph_seed in graph_seed_range:
            g = generate_graph(round(n), graph_seed, 2, 2.5, 0, 10, 10)
            for e in epsilons:
                for random_seed in random_seed_range:
                    ret = g.color(e, factoring, random_seed)
                    data.append(ret)

            df = pd.DataFrame(data=data, columns=data_labels)
            df.to_csv(name)

#              "mul" or "exp" |
#                             v
data_gen(range(3), range(3), "mul", [0.25, 0.5, 1], [10000], "datasets/new_data.csv")
