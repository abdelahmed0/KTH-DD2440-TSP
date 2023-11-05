#include <iostream>
#include <vector>
#include <limits>
#include <chrono>
#include <random>
#include <numeric>
#include <algorithm>
#include <fstream>

using namespace std;

#define TESTING 0

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
inline Matrix createNearestNeighborMatrix(Matrix& d, const size_t K_NEAREST) {
    size_t n = d.rows();
    size_t K = min(K_NEAREST, n);
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

/**
 * Performs 2-opt local optimization on a tour
 * Implementation of fast-2-opt explained in paper 
 * "Large-Step Markov Chains for the Traveling Salesman Problem"
 * Time complexity: O(n)
 * 
 * @return True if an improvement is made, false otherwise. FIXME: always true for some reason
 */
inline bool fastTwoOpt(Matrix& d, Matrix& nbhd, 
                       vector<uint32_t>& tour, vector<uint32_t> whichSlot, 
                       uint32_t minLink, uint32_t& maxLink) {
    size_t N = d.rows();
    bool improved = false;

    // We want to consider swaps of edges (n1, n2) and (m1, m2) which are currently in tour
    // for each (n1, n2)
    for (size_t n1_i = 0; n1_i < N; ++n1_i) {
        size_t n2_i = (n1_i + 1) % N;
        uint32_t n1 = tour[n1_i];
        uint32_t n2 = tour[n2_i];
        
        // for each (m1, m2) where m1 is choosen as the k-closest neighbor of n1
        for (size_t k = 0; k < nbhd.cols(); ++k) {
            uint32_t m1_i = whichSlot[nbhd.at(n1, k)];
            size_t m2_i = (m1_i + 1) % N;
            uint32_t m1 = tour[m1_i];
            uint32_t m2 = tour[m2_i];

            // if lower bound on new length is greater than upper bound of old length
            if (d.at(n1, m1) + minLink > d.at(n1, n2) + maxLink) {
                break; // go to next n1
            }
            if (d.at(n1, m1) + d.at(n2, m2) < d.at(n1, n2) + d.at(m1, m2)) {
                // make swap
                improved = true;
                reverseTourSegment(tour, whichSlot, n2_i % N, m1_i);
                maxLink = max(maxLink, max(d.at(n1, m1), d.at(n2, m2)));
                break; // go to next n1
            }
        }
    }
    return improved;
}


int main(int argc, char *argv[]) {
    auto startTime = now();
    const uint16_t timeLimit = 1900;
    const uint16_t twoOptTime = 200;
    const size_t K_NEAREST = 30;

    #if TESTING
        fstream testFile;
        testFile.open(argv[1], ios::in);
        Matrix distanceMatrix = createDistMatrixFromInput(testFile);
    #else
        Matrix distanceMatrix = createDistMatrixFromInput(cin);
    #endif

    Matrix nbhd = createNearestNeighborMatrix(distanceMatrix, K_NEAREST);
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

    bool twoOptMinima = false;
    while (!timeOver(startTime, twoOptTime)) {
        if (!fastTwoOpt(distanceMatrix, nbhd, tour, whichSlot, minLink, maxLink)) {
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
