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

static bool timeOver(chrono::steady_clock::time_point start_time, uint16_t ms_to_run) {
    return chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() >= ms_to_run;
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
 * For nbhd[i][j] = k , city k i s the jth closest city to city i
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
 * Performs fast the 2-opt local optimization algorithm from the paper
 * "Large-Step Markov Chains for the Traveling Salesman Problem"
 * O(N) checkout time
 */
inline void fastTwoOpt(Matrix& d, Matrix& nbhd,
                       vector<uint32_t>& tour, vector<uint16_t>& whichSlot,
                       uint32_t minLink, uint32_t& maxLink,
                       chrono::steady_clock::time_point startTime, uint32_t timeLimit) {
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

/**
 * Sorts three edges according to their position in the tour, given by the indices
 */
inline void sortInTourOrder(vector<uint32_t>& cities, vector<uint16_t>& ind) {
    if (ind[0] < ind[4] && ind[0] > ind[2] 
        || ind[0] > ind[2] && ind[2] > ind[4]
        || ind[0] < ind[4] && ind[2] > ind[4]) {
            vector<uint32_t> tmp_c = cities;
            cities[0] = tmp_c[2]; cities[1] = tmp_c[3];
            cities[2] = tmp_c[0]; cities[3] = tmp_c[1];
            // cities[4] and cities[5] do not change
            vector<uint16_t> tmp_i = ind;
            ind[0] = tmp_i[2]; ind[1] = tmp_i[3];
            ind[2] = tmp_i[0]; ind[3] = tmp_i[1];
            // and their indices neither
        }
}

/**
 * Performs fast the 3-opt local optimization algorithm from the paper
 * "Large-Step Markov Chains for the Traveling Salesman Problem"
 * O(N) checkout time
 */
inline void fastThreeOpt(Matrix& d, Matrix& nbhd, 
                     vector<uint32_t>& tour, vector<uint16_t>& whichSlot, 
                     uint32_t minLink, uint32_t& maxLink, 
                     chrono::steady_clock::time_point startTime, uint32_t timeLimit) {
    const uint16_t N = d.rows();

    bool improved = true;
    while (improved) {
        improved = false;

        // for each edge (m1, n1)
        for (uint16_t n1_i = 0; n1_i < N; ++n1_i) {
            uint16_t m1_i = (n1_i - 1 + N) % N;
            uint32_t m1 = tour[m1_i];
            uint32_t n1 = tour[n1_i];

            // moved here from outer loop because of timelimit issues for k > 40
            if (timeOver(startTime, timeLimit))
                return; 

            // for each edge (m2, n2) where n2 is choosen as the  k-nearest neighbor of m1
            for (uint16_t k1 = 0; k1 < nbhd.cols(); ++k1) {
                uint16_t n2_i = whichSlot[nbhd.at(m1, k1)];
                uint16_t m2_i = (n2_i + N - 1) % N; 
                uint32_t m2 = tour[m2_i];
                uint32_t n2 = tour[n2_i];

                bool newBestMove = false;

                if (d.at(m1, n2) + 2 * minLink > d.at(m1, n1) + 2 * maxLink)
                    break; // next (m1, n1) edge

                if (m1 == m2 
                    || m2 == n1
                    || d.at(m1, n2) + 2 * minLink > d.at(m1, n1) + d.at(m2, n2) + maxLink)
                    continue; // next (m2, n2) edge

                // for each (m3, n3) where n3 is chosen as the k-closest neighbor of m2
                for (uint16_t k2 = 0; k2 < nbhd.cols(); ++k2) {
                    uint16_t n3_i = whichSlot[nbhd.at(m1, k2)];
                    uint16_t m3_i = (n3_i + N - 1) % N;
                    uint32_t m3 = tour[m3_i];
                    uint32_t n3 = tour[n3_i];

                    if (d.at(m1, n2) + d.at(n1, n3) + minLink > d.at(m1, n1) + d.at(m2, n2) + maxLink)
                        break; // next (m2, n2) edge

                    if (n3 == n2 || m3 == m1
                        || n3 == m2 ||  m3 == n2
                        || m3 == n1)
                        continue; // next (m3, n3) edge

                    // sort cities according to their occurence in the tour
                    // and keep track of their indices for tour segment reversal
                    vector<uint32_t> C = {m1, n1, m2, n2, m3, n3};
                    vector<uint16_t> I = {m1_i, n1_i, m2_i, n2_i, m3_i, n3_i};
                    sortInTourOrder(C, I);

                    uint32_t oldDelta = d.at(C[0], C[1]) + d.at(C[2], C[3]) + d.at(C[4], C[5]);

                    // try the two topologically different moves
                    // which correspond to four different ordered cases 
                    // cases that repeat cities degenerate to twoOpt cases
                    if (d.at(C[3], C[0]) + d.at(C[5], C[1]) + d.at(C[2], C[4]) < oldDelta) {
                        newBestMove = true;
                        reverseTourSegment(tour, whichSlot, I[5], I[0]);
                        reverseTourSegment(tour, whichSlot, I[3], I[4]);
                        maxLink = max(maxLink, max(max(d.at(C[3], C[0]), d.at(C[5], C[1])), d.at(C[2], C[4])));
                    }  else if (d.at(C[1], C[3]) + d.at(C[4], C[0]) + d.at(C[5], C[2]) < oldDelta) {
                        newBestMove = true;
                        reverseTourSegment(tour, whichSlot, I[5], I[0]);
                        reverseTourSegment(tour, whichSlot, I[1], I[2]);
                        maxLink = max(maxLink, max(max(d.at(C[1], C[3]), d.at(C[4], C[0])), d.at(C[5], C[2])));
                    } else if (d.at(C[0], C[2]) + d.at(C[1], C[4]) + d.at(C[3], C[5]) < oldDelta) {
                        newBestMove = true;
                        reverseTourSegment(tour, whichSlot, I[1], I[2]);
                        reverseTourSegment(tour, whichSlot, I[3], I[4]);
                        maxLink = max(maxLink, max(max(d.at(C[0], C[2]), d.at(C[1], C[4])), d.at(C[3], C[5])));
                    } else if (d.at(C[1], C[4]) + d.at(C[3], C[0]) + d.at(C[5], C[2]) < oldDelta) {
                        newBestMove = true;
                        reverseTourSegment(tour, whichSlot, I[0], I[5]);
                        reverseTourSegment(tour, whichSlot, I[1], I[2]);
                        reverseTourSegment(tour, whichSlot, I[3], I[4]);
                        maxLink = max(maxLink, max(max(d.at(C[1], C[4]), d.at(C[3], C[0])), d.at(C[5], C[2])));
                    }

                    if (newBestMove) {
                        improved = true;
                        break;
                    }
                }
                if (newBestMove) {
                    break;
                }
            }
        }
    }
}

inline void setBestTour(Matrix& distanceMatrix, vector<uint32_t>& bestTour, uint32_t& bestCost, vector<uint32_t>& newTour) {
    uint32_t length = tourLength(newTour, distanceMatrix);
    if (length < bestCost) {
        bestCost = length;
        bestTour = newTour;
    }
}

/**
 * 4-opt move as explained in 
 * "An Effective Implementation of K-opt Moves for the Lin-Kernighan TSP Heuristic"
 */
inline vector<uint32_t> randomFourOptMove(const vector<uint32_t>& tour) {
    vector<uint32_t> newTour;
    const uint16_t distSize = tour.size() / 4;

    // sample random cities
    const uint16_t m1 = rand() % distSize;
    const uint16_t m2 = m1 + 1 + rand() % distSize;
    const uint16_t m3 = m2 + 1 + rand() % distSize;

    // create new tour with double bridging old tour
    newTour.insert(newTour.end(), tour.begin() + m3, tour.end());
    newTour.insert(newTour.end(), tour.begin() + m1, tour.begin() + m2);
    newTour.insert(newTour.end(), tour.begin() + m2, tour.begin() + m3);
    newTour.insert(newTour.end(), tour.begin(), tour.begin() + m1);

    return newTour;
}

int main(int argc, char *argv[]) {
    auto startTime = chrono::steady_clock::now();
    const uint16_t timePerTwoOpt = 10;
    const uint16_t timePerIteration = 30;
    const uint16_t timeLimit = 2000 - timePerIteration;
    
    // k=20 is a typical value according to "The Traveling Salesman Problem: A Case Study in Local Optimization"
    // for random instances, k=15 seems to be optimal
    const uint16_t K_NEAREST = 20;

#if TESTING
    fstream testFile;
    testFile.open(argv[1], ios::in);
    Matrix distanceMatrix = createDistMatrixFromInput(testFile);
#else
    Matrix distanceMatrix = createDistMatrixFromInput(cin);
#endif

    if (distanceMatrix.rows() == 1) {
        std::cout << 0 << std::endl;
        return 0;
    }

    Matrix nbhd = createNearestNeighborMatrix(distanceMatrix, K_NEAREST);
    uint16_t n = distanceMatrix.rows();

    vector<uint32_t> tour = greedy(distanceMatrix);

#if TESTING
    cout << "Greedy length: " << tourLength(tour, distanceMatrix) << endl;
#endif

    vector<uint32_t> bestTour = tour;
    uint32_t bestCost = tourLength(tour, distanceMatrix);

    while (!timeOver(startTime, timeLimit)){
        auto iterationStartTime = chrono::steady_clock::now();
        // Vector to keep track of where cities are in the current tour
        // Used for fast-2-opt and fast-3-opt
        vector<uint16_t> whichSlot(n);
        for (uint16_t i = 0; i < n; ++i) {
            whichSlot[tour[i]] = i;
        }
        const uint32_t minLink = distanceMatrix.min();
        uint32_t maxLink = *max_element(begin(tour), end(tour));

        fastTwoOpt(distanceMatrix, nbhd, tour, whichSlot, minLink, maxLink, iterationStartTime, timePerTwoOpt);
        setBestTour(distanceMatrix, bestTour, bestCost, tour);

        fastThreeOpt(distanceMatrix, nbhd, tour, whichSlot, minLink, maxLink, iterationStartTime, timePerIteration);
        setBestTour(distanceMatrix, bestTour, bestCost, tour);
        
        // random tour perturbation by using a 4-Opt move
        tour = randomFourOptMove(tour);
    }


#if TESTING
    cout << "Tourlength:    " << bestCost;
#else
    for (int i = 0; i < bestTour.size(); ++i) {
        cout << bestTour[i] << endl;
    }
#endif    

    return 0;
}
