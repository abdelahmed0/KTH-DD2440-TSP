#include <iostream>
#include <vector>
#include <cmath>

#define MAX_N 1000

typedef double num_t;

typedef struct {
    num_t x;
    num_t y;
} vec_t;


int n;
vec_t points[MAX_N] = {0};

num_t distance(vec_t& a, vec_t& b) {
    num_t x = a.x - b.x;
    num_t y = a.y - b.y;
    return std::round(std::sqrt(x * x + y * y));
}

void greedy_tour();

int main() {
    std::cin >> n;

    for (int i = 0; i < n; ++i) {
        num_t x, y;
        std::cin >> x >> y;
        points[i] = {x, y};
    }

    greedy_tour();

    return 0;
}

void greedy_tour() {
    int tour[n];
    bool used[n];

    for (int i = 1; i < n; ++i) {
        tour[i] = i;
        used[i] = false;
    }

    tour[0] = 0;
    used[0] = true;

    for (int i = 1; i < n; ++i) {
        int best = -1;
        for (int j = 0; j < n; ++j) {
            vec_t tour_v = points[tour[i - 1]];
            vec_t j_v = points[j];
            vec_t best_v = points[best];
            if ((!used[j]) && (best == -1 || distance(tour_v, j_v) < distance(tour_v, best_v))) {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }

    for (int i = 0; i < n; ++i) {
        std::cout << tour[i] << std::endl;
    }
}