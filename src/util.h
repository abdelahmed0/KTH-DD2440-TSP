#ifndef TSP_UTIL_H
#define TSP_UTIL_H

#include <cstdint>
#include <vector>
#include <array>
#include <set>

#define NEAREST 20

typedef double num_t;
typedef uint32_t length_t;

typedef struct {
    int first;
    int second;
} edge_t;

bool operator<(const edge_t& lhs, const edge_t& rhs);
bool operator==(const edge_t& lhs, const edge_t& rhs);

edge_t edge(int from, int to);

bool set_contains(const std::set<edge_t>& set, edge_t edge);


class Matrix {
public:
    explicit Matrix(int n);
    int dim() const;
    length_t& at(int i, int j);

private:
    int n;
    std::vector<length_t> entries;
};

class Neighbours {
public:
    explicit Neighbours(int n);
    void find_neighbours(int idx, Matrix &distances);
    std::vector<int>& at(int idx);

private:
    std::vector<std::vector<int>> entries;
};


class Tour {
public:
    explicit Tour(int n);
    length_t length(Matrix& distances);
    bool is_valid();

    int& operator[](int i);
    // returns the index of the node
    int index_of(int node);
    // returns predecessor node of i
    int predecessor(int i);
    // returns the successor node of i
    int successor(int i);
    // returns both the predecessor and successor node of i
    std::array<int, 2> adjacent(int i);
    // returns true iff there is an (undirected) edge between t1 and t2
    bool are_connected(int t1, int t2);

private:
    int n;
    std::vector<int> tour;
};

#endif //TSP_UTIL_H
