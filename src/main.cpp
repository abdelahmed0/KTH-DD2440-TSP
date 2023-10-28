#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <set>
#include "util.h"
#include "lin_kernighan.h"

#define RUNTIME 1900

static inline std::chrono::time_point<std::chrono::high_resolution_clock> now() {
    return std::chrono::high_resolution_clock::now();
}

static bool is_time_over(std::chrono::high_resolution_clock::time_point start_time, uint16_t ms_to_run) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(now() - start_time).count() >= ms_to_run;
}

void greedy_tour(Matrix& distances, Tour& tour);

int main() {
    auto start = now();

    int n;
    std::cin >> n;

    std::vector<std::pair<num_t, num_t>> points;
    for (int i = 0; i < n; ++i) {
        num_t x, y;
        std::cin >> x >> y;
        points.emplace_back(x, y);
    }

    // populate distances
    Matrix distances(n);
    for (int i = 0; i < n; ++i) {
        auto& a = points[i];
        for (int j = 0; j < n; ++j) {
            auto& b = points[j];
            num_t x = a.first - b.first;
            num_t y = a.second - b.second;
            // NOTE possibly change to other cast
            distances.at(i, j) = static_cast<length_t>(std::round(std::sqrt(x * x + y * y)));
        }
    }

    // compute neighborhood for points
    Neighbours neighbours(n);
    for (int i = 0; i < n; ++i) {
        neighbours.find_neighbours(i, distances);
    }

    Tour tour(n);
    greedy_tour(distances, tour);

    LK lk(tour, distances, neighbours);
    while (!is_time_over(start, RUNTIME)) {
        if (!lk.step()) {
            break;
        }
    }

    for (int i = 0; i < n; ++i) {
        std::cout << lk.get_tour()[i] << std::endl;
    }

    return 0;
}

void greedy_tour(Matrix& distances, Tour& tour) {
    int n = distances.dim();
    bool used[n];

    for (int i = 1; i < n; ++i) {
        used[i] = false;
    }

    tour[0] = 0;
    used[0] = true;

    for (int i = 1; i < n; ++i) {
        int best = -1;
        for (int j = 0; j < n; ++j) {
            if ((!used[j]) && (best == -1 || distances.at(tour[i - 1], j) < distances.at(tour[i - 1], best))) {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }
}