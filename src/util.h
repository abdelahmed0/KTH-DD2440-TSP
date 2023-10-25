#ifndef TSP_UTIL_H
#define TSP_UTIL_H

// number of neighbours stored for each city (for optimization)
#define NEAREST 20

typedef double num_t; // in case we want to change to float's later for (possibly) faster runtime

// data structure to store 2d coordinates
typedef struct {
    num_t x;
    num_t y;
} vec_t;

typedef struct city_t {
    vec_t pos{};
    city_t* nearest[NEAREST] = {nullptr};
} city_t;

num_t distance(vec_t& a, vec_t& b);
// finds the NEAREST number of neighbours and updates city->nearest with them
void compute_neighbours(int n, city_t world[], city_t& city);

#endif //TSP_UTIL_H

