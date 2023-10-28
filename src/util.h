#ifndef TSP_UTIL_H
#define TSP_UTIL_H

#include <cstdint>
#include <vector>
#include <array>

#define NEAREST 20

typedef double num_t;
typedef uint32_t length_t;
typedef std::pair<int, int> edge_t;

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

    Tour update(std::set<edge_t>& added, std::set<edge_t>& removed);

private:
    int n;
    std::vector<int> tour;
};

#endif //TSP_UTIL_H
