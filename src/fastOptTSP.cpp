#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <numeric>
#include <algorithm>
#include <fstream>

using namespace std;

static bool timeOver(const chrono::steady_clock::time_point start_time, const chrono::milliseconds end_time) {
    chrono::milliseconds elapsed_time = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_time);
    return elapsed_time >= end_time;
}

class Matrix {
public:
    Matrix(uint16_t N, uint16_t M) : n(N), m(M), data(N * M) {}
    inline uint16_t rows() { return n; }
    inline uint16_t cols() { return m; }
    uint32_t& at(uint16_t i, uint16_t j) { return data[i * m + j]; }
    vector<uint32_t>::iterator rowIterator(uint16_t i) { return data.begin() + i * m; }
    uint32_t min() {return *min_element(begin(data), end(data));}

private:
    uint16_t n;
    uint16_t m;
    vector<uint32_t> data;
};

inline Matrix createDistMatrixFromInput(istream& in) {
    uint16_t n;
    in >> n;
    Matrix distanceMatrix(n, n);
    
    vector<double> x(n);
    vector<double> y(n);
    for (int line = 0; line < n; ++line) {
        in >> x[line] >> y[line];
    }

    for (uint16_t r = 0; r < n; ++r) {
        for (uint16_t c = 0; c < n; ++c) {
            distanceMatrix.at(r, c) = round(sqrt(pow(x[r] - x[c], 2) + pow(y[r] - y[c], 2)));
        }
    }

    return distanceMatrix;
}

/**
 * Calculate K-nearest neighbor matrix
 * Where nbhd[i][k] is the k-nearest city to i
 */
inline Matrix createNearestNeighborMatrix(Matrix& d, const uint16_t K_NEAREST) {
    uint16_t n = d.rows();
    uint16_t K = min(K_NEAREST, n);
    Matrix nbhd(n, K);
    vector<uint32_t> nbhdRow(n-1);
        
    for (uint16_t i = 0; i < n; ++i) {
        iota(nbhdRow.begin(), nbhdRow.begin() + i, 0);
        iota(nbhdRow.begin() + i, nbhdRow.end(), i+1);

        stable_sort(nbhdRow.begin(), nbhdRow.end(),
            [&](uint32_t j, uint32_t k) {
                return d.at(i, j) < d.at(i, k);
            }
        );
        for (uint16_t k = 0; k < K; ++k) {
            nbhd.at(i, k) = nbhdRow[k];
        }
    }
    
    return nbhd;
}

inline void reverseTourSegment(vector<uint32_t>& tour, vector<uint16_t>& whichSlot, uint16_t start, uint16_t end) {
    while (start < end) {
        swap(tour[start], tour[end]);
        whichSlot[tour[start]] = start;
        whichSlot[tour[end]] = end;
        start++;
        end--;
    }
}

inline uint32_t tourLength(vector<uint32_t>& tour, Matrix dist) { 
    uint32_t length = 0;
    uint16_t n = tour.size();
    for (uint16_t i = 0; i < n; ++i) {
        length += dist.at(tour[i], tour[(i+1) % n]);
    }
    return length;
}

inline vector<uint32_t> greedy(Matrix& m) {
    uint16_t n = m.rows();
    vector<uint32_t> tour(n);
    vector<bool> used(n, false);
    tour[0] = 0;
    used[0] = true;

    for (uint16_t i = 1; i < n; ++i) {
        uint16_t currentCity = tour[i - 1];
        uint16_t nearestNeighbor = -1;
        uint32_t minDistance = -1;

        for (uint16_t j = 0; j < n; ++j) {
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
 * Performs fast the 2-opt local optimization algorithm from the paper
 * "Large-Step Markov Chains for the Traveling Salesman Problem"
 * O(N) checkout time
 */
inline void fastTwoOpt(Matrix& d, Matrix& nbhd,
                       vector<uint32_t>& tour, vector<uint16_t>& whichSlot,
                       uint32_t minLink, uint32_t& maxLink,
                       chrono::steady_clock::time_point startTime, chrono::milliseconds timeLimit) {
    uint16_t N = d.rows();
    bool improved = true;

    while (improved) {
        if (timeOver(startTime, timeLimit)) break;

        improved = false;
        // We want to consider swaps of edges (m1, n1) and (m2, n2) which are currently in tour
        // for each (m1, n1)
        for (uint16_t n1_i = 0; n1_i < N; ++n1_i) {
            uint16_t m1_i = (n1_i - 1 + N) % N;
            uint32_t n1 = tour[n1_i];
            uint32_t m1 = tour[m1_i];

            // for each (m2, n2) where m2 is choosen as the k-closest neighbor of m1
            for (uint16_t k = 0; k < nbhd.cols(); ++k) {
                uint32_t n2_i = whichSlot[nbhd.at(n1, k)];
                uint16_t m2_i = (n2_i - 1 + N) % N;
                uint32_t n2 = tour[n2_i];
                uint32_t m2 = tour[m2_i];

                // if lower bound on new length is greater than upper bound of old length
                if (d.at(n1, n2) + minLink > d.at(m1, n1) + maxLink) {
                    break; // go to next m1
                }
                // try the move
                if (d.at(m1, m2) + d.at(n1, n2) < d.at(m1, n1) + d.at(m2, n2)) {
                    improved = true;
                    reverseTourSegment(tour, whichSlot, n1_i, m2_i);
                    maxLink = max(maxLink, max(d.at(m1, m2), d.at(n1, n2)));
                    break; // go to next m1
                }
            }
        }
    }
}

inline bool setBestTour(Matrix& distanceMatrix, vector<uint32_t>& bestTour, uint32_t& bestCost, vector<uint32_t>& newTour) {
    uint32_t length = tourLength(newTour, distanceMatrix);
    if (length < bestCost) {
        bestCost = length;
        bestTour = newTour;
        return true;
    }
    return false;
}

/**
 * Large-step markov-chain tour perturbation
 * 4-opt move as explained in 
 * "An Effective Implementation of K-opt Moves for the Lin-Kernighan TSP Heuristic"
 * TODO: try (fast) pseudo random numbers with seed
 */
inline vector<uint32_t> randomFourOptMove(const vector<uint32_t>& tour) {
    vector<uint32_t> newTour;
    const uint16_t distSize = tour.size() / 3;

    // sampling three random cities gives us 4 toursegments (start and end included)
    // [start, m1] [m1, m2] [m2, m3] [m3, end]
    const uint16_t m1 = rand() % distSize;
    const uint16_t m2 = m1 + rand() % distSize;
    const uint16_t m3 = m2 + rand() % distSize;

    // create new tour with double bridging old tour
    // which corresponds to three link exchanges and results in tour
    // [start, m1] [m3, end] [m1, m2] [m2, m3]  
    newTour.insert(newTour.end(), tour.begin() + m3, tour.end());
    newTour.insert(newTour.end(), tour.begin() + m1, tour.begin() + m2);
    newTour.insert(newTour.end(), tour.begin() + m2, tour.begin() + m3);
    newTour.insert(newTour.end(), tour.begin(), tour.begin() + m1);

    return newTour;
}

int main() {
    const auto startTime = chrono::steady_clock::now();

    // k=20 is a typical value according to "The Traveling Salesman Problem: A Case Study in Local Optimization"
    // for random instances, k=15 seems to be optimal
    const uint16_t K_NEAREST = 20;
    const auto timePerTwoOpt = chrono::milliseconds(5);
    // const auto timePerIteration = chrono::milliseconds(10);
    const auto timeLimit = chrono::milliseconds(1998) - timePerTwoOpt;  
    // const auto timeLimit = chrono::milliseconds(1998) - timePerIteration;

    Matrix distanceMatrix = createDistMatrixFromInput(cin);

    if (distanceMatrix.rows() == 1) {
        std::cout << 0 << std::endl;
        return 0;
    }

    Matrix nbhd = createNearestNeighborMatrix(distanceMatrix, K_NEAREST);
    const uint16_t n = distanceMatrix.rows();
    const uint32_t minLink = distanceMatrix.min();

    vector<uint32_t> trialTour = greedy(distanceMatrix);
    uint32_t bestCost = tourLength(trialTour, distanceMatrix);

    vector<uint32_t> bestTour = trialTour;

    while (!timeOver(startTime, timeLimit)) {
        auto iterationStartTime = chrono::steady_clock::now();
        // Vector to keep track of where cities are in the current tour
        // Used for fast-opt
        vector<uint16_t> whichSlot(n);
        for (uint16_t i = 0; i < n; ++i) {
            whichSlot[trialTour[i]] = i;
        }
        uint32_t maxLink = *max_element(begin(trialTour), end(trialTour));

        bool improved = false;
        fastTwoOpt(distanceMatrix, nbhd, trialTour, whichSlot, minLink, maxLink, iterationStartTime, timePerTwoOpt);
        improved = setBestTour(distanceMatrix, bestTour, bestCost, trialTour);

        if (improved) {
            trialTour = randomFourOptMove(bestTour);
        } else if (rand() % 2 > 0) {
            // Tried accept according to Boltzmann distribution but had a big negative impact on performance
            trialTour = randomFourOptMove(bestTour);
        } else {
            trialTour = randomFourOptMove(trialTour);
        }
    }

    for (int i = 0; i < bestTour.size(); ++i) {
        cout << bestTour[i] << endl;
    }

    return 0;
}
