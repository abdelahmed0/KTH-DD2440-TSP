#include "graph.h"

void add_edge(graph_t& graph, node_t *from, node_t *to) {
    if(from->successor != nullptr) {
        remove_edge(graph, from, from->successor);
    }

    if(to->predecessor != nullptr) {
        remove_edge(graph, to->predecessor, to);
    }

    from->successor = to;
    to->predecessor = from;

    graph.length += distance(from->value.pos, to->value.pos);
}

void remove_edge(graph_t& graph, node_t *from, node_t *to) {
    if(from->successor == to && from == to->predecessor) {
        from->successor = nullptr;
        to->predecessor = nullptr;

        graph.length -= distance(from->value.pos, to->value.pos);
    }
}

void from_tour(graph_t &graph, const int *tour) {
    for (int i = 0; i < graph.n; ++i) {
        int j = (i + 1) % graph.n;
        int curr = tour[i];
        int next = tour[j];
        add_edge(graph, &graph.nodes[curr], &graph.nodes[next]);
    }
}

void to_tour(graph_t &graph, int *tour) {
    node_t* node = &graph.nodes[0];
    for (int i = 0; i < graph.n; ++i) {
        tour[i] = node->idx;
        node = node->successor;
    }
}
