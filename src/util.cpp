#include <set>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "util.h"

Matrix::Matrix(int n) : n(n), entries(n * n) {}

int Matrix::dim() const {
    return n;
}

length_t& Matrix::at(int i, int j) {
    return entries[i * n + j];
}

std::vector<int> &Neighbours::at(int idx) {
    return entries[idx];
}

void Neighbours::find_neighbours(int idx, Matrix &distances) {
    std::set<std::pair<length_t, int>> ordered;

    for (int i = 0; i < distances.dim(); ++i) {
        if (i == idx) {
            continue;
        }
        ordered.emplace(distances.at(idx, i), i);
    }

    entries[idx] = std::vector<int>();
    int i = 0;
    for (const auto& [len, j] : ordered) {
        if (i >= NEAREST) {
            break;
        }
        entries[idx].push_back(j);
        i++;
    }

}

Neighbours::Neighbours(int n) : entries(n) {

}

Tour::Tour(int n) : n(n){
    for (int i = 0; i < n; ++i) {
        tour.push_back(i);
    }
}

int &Tour::operator[](int i) {
    return tour[i];
}

length_t Tour::length(Matrix& distances) {
    length_t len = 0;
    for (int i = 0; i < distances.dim(); ++i) {
        int j = (i + 1) % distances.dim();
        len += distances.at(i, j);
    }
    return len;
}

bool Tour::is_valid() {
    std::set<int> visited;
    for (int i = 0; i < n; ++i) {
        visited.insert(tour[i]);
    }
    return visited.size() == n;
}

int Tour::predecessor(int i) {
    return tour[(i - 1) % n];
}

int Tour::successor(int i) {
    return tour[(i + 1) % n];
}

std::array<int, 2> Tour::adjacent(int i) {
    return {predecessor(i), successor(i)};
}

Tour Tour::update(std::set<edge_t>& added, std::set<edge_t>& removed) {
    std::set<edge_t> edges;
    for (int i = 0; i < n; ++i) {
        edge_t e = edge(tour[i], tour[(i + 1) % n]);
        if (!set_contains(removed, e)) {
            edges.insert(e);
        }
    }
    for (auto e : added) {
        edges.insert(e);
    }

    assert(edges.size() == n);

    Tour t2(n);
    int pointer = 0;
    while (!edges.empty()) {
        edge_t e = *edges.begin();
    }

    int idx = 0;
    t2[0] = idx;

    // TODO check if works

    for (int i = 1; i < n; ++i) {
        edge_t found_edge(-1, -1);
        for (auto& e : edges) {
            if (e.first == idx) {
                idx = e.second;
                found_edge = e;
            } else if (e.second == idx) {
                idx = e.first;
                found_edge = e;
            }
        }
        t2[i] = idx;
        if (found_edge != std::pair(-1, -1)) {
            edges.erase(found_edge);
        }
    }

    return t2;
}

int Tour::index_of(int node) {
    auto it = std::find(tour.begin(), tour.end(), node);
    if (it == tour.end()) {
        std::cerr << "Accessed non-present node: " << node << std::endl;
        return -1;
    }
    return (int) (it - tour.begin());
}


edge_t edge(int from, int to) {
    if (from < to)
        return {from, to};
    return {to, from};
}

bool set_contains(const std::set<edge_t> &set, edge_t edge) {
    return set.count(edge) > 0;
}
