#ifndef TSP_LIN_KERNIGHAN_H
#define TSP_LIN_KERNIGHAN_H

#include "util.h"

class LK {
public:
    explicit LK(Tour &tour, Matrix &distances, Neighbours &neighbours);
    Tour& get_tour();
    bool step();
    bool naive();

private:
    bool generate_new_tour(Tour& t2, std::set<edge_t>& added, std::set<edge_t>& removed);
    bool recurse(Tour& T, std::set<edge_t> Xi, std::set<edge_t> Yi, int t1, int last, length_t Gi);
    int choose_yi(int last, int& start);

protected:
    Tour& tour;
    Matrix& distances;
    Neighbours& neighbours;

private:
    int iteration = 0;
    Tour out;
};

#endif //TSP_LIN_KERNIGHAN_H
