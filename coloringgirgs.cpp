#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "main.cpp"

namespace py = pybind11;

Graph generate_graph_interface(int n, int graph_seed, int d = 2, double ple = 2.5, double alpha = std::numeric_limits<double>::infinity(), double deg = 10, int threads = 10)
{
    return generate_graph(n, graph_seed, d, ple, alpha == 0 ? std::numeric_limits<double>::infinity() : alpha, deg, graph_seed, graph_seed, graph_seed, 10);
}

PYBIND11_MODULE(coloringgirgs, cg)
{
    cg.doc() = "Testing coloringgirgs";
    py::class_<Graph>(cg, "Graph")
        .def(py::init<std::vector<std::pair<int, int>>, int, int, double, int>())
        .def(py::pickle(
            [](Graph &g)
            { return py::make_tuple(g.get_edges(), g.get_n(), g.get_average_degree(), g.get_alpha(), g.get_graph_seed()); },
            [](py::tuple t)
            { return Graph(t[0].cast<std::vector<std::pair<int, int>>>(), t[1].cast<int>(), t[2].cast<int>(), t[3].cast<double>(), t[4].cast<int>()); }))
        .def("color", &Graph::color, py::call_guard<py::gil_scoped_release>());
    cg.def("generate_graph", &generate_graph_interface, py::call_guard<py::gil_scoped_release>(), "generate a new graph for calculation with the given parameters");
}
