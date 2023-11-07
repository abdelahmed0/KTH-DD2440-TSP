#ifndef TSP_LIN_KERNIGHAN_H
#define TSP_LIN_KERNIGHAN_H

#include <set>
#include "util.h"

class LK {
public:
    LK(Tour& base, Matrix &distances, Neighbours &neighbours);
    Tour& get_tour();
    bool naive(stopwatch_t& stopwatch);


private:
    bool new_tour(std::vector<edge_t>& X, std::vector<edge_t>& Y, edge_t final, int i);

protected:
    Matrix& distances;
    Neighbours& neighbours;

private:
    Tour tour1, tour2;
    int iteration = 0;

    std::set<edge_t> edges;
    std::set<int> visited;

    std::vector<edge_t> X;
    std::vector<edge_t> Y;
    std::vector<length_t> Gi;
    std::vector<int> nodes;
};

#endif //TSP_LIN_KERNIGHAN_H
