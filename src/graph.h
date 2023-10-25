#ifndef TSP_GRAPH_H
#define TSP_GRAPH_H
#include "util.h"
#include <vector>

typedef struct node_t {
    int idx = -1;
    city_t value{};
    node_t* predecessor = nullptr;
    node_t* successor = nullptr;
} node_t;

typedef struct graph_t {
    int n;
    num_t length = 0.0;
    std::vector<node_t> nodes;
} graph_t;

// adds the edge (and removes dangling edges) updates the length accordingly
void add_edge(graph_t& graph, node_t* from, node_t* to);
// removes the edge (if present) updates the length accordingly
void remove_edge(graph_t& graph, node_t* from, node_t* to);
// generates the edges according to some tour (this does not populate the nodes nor set n)
void from_tour(graph_t& graph, const int tour[]);
void to_tour(graph_t& graph, int tour[]);

#endif //TSP_GRAPH_H