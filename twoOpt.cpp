#include <iostream>
#include <vector>
#include <limits>
#include <chrono>
#include <random>
#include <numeric>
#include <algorithm>
#include <fstream>

using namespace std;

#define K_NEAREST 10
#define TESTING 1

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
    vector<uint32_t>::iterator rowIterator(size_t i) { return data.begin() + i * m; }
    uint32_t min() {return *min_element(begin(data), end(data));}

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

/**
 * Calculate K-nearest neighbor matrix
 * For nbhd[i][j] = k , city k i s the jth closest city to city i
 */
inline Matrix createNearestNeighborMatrix(Matrix& d) {
    size_t n = d.rows();
    size_t K = min((size_t)K_NEAREST, n);
    Matrix nbhd(n, K);
    vector<uint32_t> nbhdRow(n-1);
        
    for (size_t i = 0; i < n; ++i) {
        iota(nbhdRow.begin(), nbhdRow.begin() + i, 0);
        iota(nbhdRow.begin() + i, nbhdRow.end(), i+1);

        stable_sort(nbhdRow.begin(), nbhdRow.end(),
            [&](uint32_t j, uint32_t k) {
                return d.at(i, j) < d.at(i, k);
            }
        );
        for (size_t k = 0; k < K; ++k) {
            nbhd.at(i, k) = nbhdRow[k];
        }
    }
    
    return nbhd;
}

Matrix createNeighborsMatrix(Matrix& d) {
    size_t N = d.rows();
    size_t M = d.cols() - 1;
    size_t K = min(M, (size_t)K_NEAREST);
    Matrix neighbor(N, K);
    vector<uint32_t> row(M); // For sorting.

    for (size_t i = 0; i < N; ++i) {
        // Fill row with 0, 1, ..., i - 1, i + 1, ..., M - 1.
        uint16_t k = 0;
        for (size_t j = 0; j < M; ++j, ++k) {
            row[j] = (i == j) ? ++k : k;
        }
        // Sort K first elements in row by distance to i.
        partial_sort(row.begin(), row.begin() + K, row.end(),
            [&](uint16_t j, uint16_t k) {
                return d.at(i, j) < d.at(i, k);
            }
        );
        // Copy first K elements (now sorted) to neighbor matrix.
        for (size_t k = 0; k < K; ++k) {
            neighbor.at(i, k) = row[k];
        }
    }
    return neighbor;
}

/**
 * Requires start to be smaller than end
 */
inline void reverseTourSegment(vector<uint32_t>& tour, vector<uint32_t>& whichSlot, size_t start, size_t end) {
    while (start < end) {
        swap(tour[start], tour[end]);
        whichSlot[tour[start]] = start;
        whichSlot[tour[end]] = end;
        start++;
        end--;
    }
}

inline uint32_t tourLength(vector<uint32_t>& tour, Matrix dist) { // TODO: check if whichSlot is needed
    uint32_t length = 0;
    size_t n = tour.size();
    for (size_t i = 0; i < n; ++i) {
        length += dist.at(tour[i], tour[(i+1) % n]);
    }
    return length;
}

inline vector<uint32_t> greedy(Matrix& m) {
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

inline bool twoOpt(Matrix& d,
        Matrix& neighbor, vector<uint32_t>& tour, vector<uint32_t> &position,
        uint32_t min, uint32_t& maxLink) {
    size_t N = d.rows(); // Number of cities.

    // Candidate edges uv, wz and their positions in tour.
    uint16_t u, v, w, z;
    size_t u_i, v_i, w_i, z_i;

    bool improved = false;
    
    // For each edge uv.
    for (u_i = 0, v_i = 1; u_i < N; ++u_i, ++v_i) {
        u = tour[u_i];
        v = tour[v_i % N];

        // For each edge wz (w k:th closest neighbor of u).
        for (size_t k = 0; k < neighbor.cols(); ++k) {
            w_i = position[neighbor.at(u, k)];
            z_i = w_i + 1;
            w = tour[w_i];
            z = tour[z_i % N];

            if (v == w || z == u) {
                continue; // Skip adjacent edges.
            }

            // d[u][w] + min is a lower bound on new length.
            // d[u][v] + max is an upper bound on old length.
            if (d.at(u, w) + min > d.at(u, v) + maxLink) {
                break; // Go to next edge uv.
            }

            if (d.at(u, w) + d.at(v, z) < d.at(u, v) + d.at(w, z)) {
                //   --u w--        --u-w->
                //      X     ===>
                //   <-z v->        <-z-v--
                // (vector<uint32_t>& tour, vector<uint32_t>& whichSlot, size_t start, size_t end) 
                reverseTourSegment(tour, position, v_i % N, w_i);
                maxLink = max(maxLink, max(d.at(u, w), d.at(v, z)));
                improved = true;
                break;
            }
        }
    }
    return improved;
}

int main(int argc, char *argv[]) {
    auto startTime = now();
    const uint16_t timeLimit = 200;

    #if TESTING
        fstream testFile;
        testFile.open(argv[1], ios::in);
        Matrix distanceMatrix = createDistMatrixFromInput(testFile);
    #else
        Matrix distanceMatrix = createDistMatrixFromInput(cin);
    #endif

    Matrix nbhd = createNearestNeighborMatrix(distanceMatrix);
    size_t n = distanceMatrix.rows();

    vector<uint32_t> tour = greedy(distanceMatrix);

    // Vector to keep track of where cities are in the current tour
    // Used for fast-2-opt and fast-3-opt
    vector<uint32_t> whichSlot(n);
    for (uint32_t i = 0; i < n; ++i) {
        whichSlot[tour[i]] = i;
    }
    const uint32_t minLink = distanceMatrix.min();
    uint32_t maxLink = *max_element(begin(tour), end(tour));

    // bool twoOptMinima = false;
    while (!timeOver(startTime, timeLimit)) {
        if (!twoOpt(distanceMatrix, nbhd, tour, whichSlot, minLink, maxLink)) {
            // twoOptMinima = true; // uncomment to activate 3-opt after 2-opt
            break;
        }
    }

    // if (twoOptMinima) { //TODO: threeOpt implementation
    //     while (true) {
    //         if (!threeOptTimedTSP(distanceMatrix, tour, startTime, timeLimit)) {
    //             break;
    //         }
    //     }
    // }

    for (int i = 0; i < tour.size(); ++i) {
        cout << tour[i] << endl;
    }
    #if TESTING
        cout << "Tourlength: " << tourLength(tour, distanceMatrix) << endl;
    #endif

    return 0;
}
