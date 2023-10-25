#include <iostream>
#include <vector>
#include <limits>
#include <chrono>
#include <random>

using namespace std;

static inline chrono::time_point<chrono::high_resolution_clock> now() {
    return chrono::high_resolution_clock::now();
}

static bool timeOver(chrono::high_resolution_clock::time_point start_time, uint16_t ms_to_run) {
    return chrono::duration_cast<chrono::milliseconds>(now() - start_time).count() >= ms_to_run;
}

class Matrix {
public:
    Matrix(size_t N, size_t M) : n(N), m(M), data(N * M) {}
    inline size_t rows() { return n; }
    inline size_t cols() { return m; }
    uint32_t& at(size_t i, size_t j) { return data[i * m + j]; }

private:
    size_t n;
    size_t m;
    vector<uint32_t> data;
};

inline Matrix createDistMatrixFromInput(istream& in) {
    size_t n;
    in >> n;
    Matrix distanceMatrix(n, n);
    
    vector<double> x(n);
    vector<double> y(n);
    for (int line = 0; line < n; ++line) {
        in >> x[line] >> y[line];
    }

    for (size_t r = 0; r < n; ++r) {
        for (size_t c = 0; c < n; ++c) {
            distanceMatrix.at(r, c) = round(sqrt(pow(x[r] - x[c], 2) + pow(y[r] - y[c], 2)));
        }
    }

    return distanceMatrix;
}

inline void reverseTourSegment(vector<uint32_t>& tour, size_t start, size_t end) {
    while (start < end) {
        swap(tour[start], tour[end]);
        start++;
        end--;
    }
}

inline uint32_t tourLength(vector<uint32_t>& tour, Matrix dist) {
    uint32_t length = 0;
    size_t n = tour.size();
    for (size_t i = 0; i < n; ++i) {
        length += dist.at(tour[i], tour[(i+1) % n]);
    }
    return length;
}

inline vector<uint32_t> greedyTSP(Matrix& m) {
    int n = m.rows();
    vector<uint32_t> tour(n);
    vector<bool> used(n, false);
    tour[0] = 0;
    used[0] = true;

    for (int i = 1; i < n; ++i) {
        int currentCity = tour[i - 1];
        int nearestNeighbor = -1;
        uint32_t minDistance = numeric_limits<uint32_t>::max();

        for (int j = 0; j < n; ++j) {
            if (!used[j] && m.at(currentCity, j) < minDistance) {
                minDistance = m.at(currentCity, j);
                nearestNeighbor = j;
            }
        }

        tour[i] = nearestNeighbor;
        used[nearestNeighbor] = true;
    }

    return tour;
}

inline bool twoOptTSP(Matrix& d, vector<uint32_t>& tour) {
    int n = d.rows();
    
    // 2-Opt local optimization
    bool improved = false;
    for (int i = 0; i < n - 2; i++) {
        for (int j = i + 2; j < n; j++) {
            int cityA = tour[i];
            int cityB = tour[i + 1];
            int cityC = tour[j];
            int cityD = tour[(j + 1) % n];

            uint32_t old_distance = d.at(cityA, cityB) + d.at(cityC, cityD);
            uint32_t new_distance = d.at(cityA, cityC) + d.at(cityB, cityD);

            if (new_distance < old_distance) {
                // Perform the 2-opt swap
                reverseTourSegment(tour, i + 1, j);
                improved = true;
            }
        }
    }
    return improved;
}

int main() {
    auto startTime = now();

    Matrix distanceMatrix = createDistMatrixFromInput(cin);
    vector<uint32_t> tour = greedyTSP(distanceMatrix);

    while (!timeOver(startTime, 1980) && twoOptTSP(distanceMatrix, tour)) {
        continue;
    }

    for (int i = 0; i < tour.size(); ++i) {
        std::cout << tour[i] << std::endl;
    }

    return 0;
}
