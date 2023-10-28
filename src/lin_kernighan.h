#ifndef TSP_LIN_KERNIGHAN_H
#define TSP_LIN_KERNIGHAN_H

#include "util.h"

class LK {
public:
    explicit LK(Tour &tour, Matrix &distances, Neighbours &neighbours);
    bool step();
    Tour& get_tour();

private:
    bool chooseX(int t1, int t2, length_t& gain, std::set<edge_t>& removed_edges, std::set<edge_t>& added_edges);
    bool chooseY(int t1, int t2, length_t& gain, std::set<edge_t>& removed_edges, std::set<edge_t>& added_edges);

protected:
    Tour& tour;
    Matrix& distances;
    Neighbours& neighbours;

private:
    int iteration = 0;
    Tour out;
};

#endif //TSP_LIN_KERNIGHAN_H
