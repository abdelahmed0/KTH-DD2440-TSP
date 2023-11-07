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
 * Performs fast the 2-opt local optimization algorithm from the paper
 * "Large-Step Markov Chains for the Traveling Salesman Problem"
 * Time complexity: O(n)
 */
inline void fastTwoOpt(Matrix& d, Matrix& nbhd,
                       vector<uint32_t>& tour, vector<uint32_t>& whichSlot,
                       uint32_t minLink, uint32_t& maxLink,
                       chrono::steady_clock::time_point startTime, uint32_t timeLimit) {
    size_t N = d.rows();
    bool improved = true;

    while (improved) {
        if (timeOver(startTime, timeLimit)) break;

        improved = false;
        // We want to consider swaps of edges (m1, n1) and (m2, n2) which are currently in tour
        // for each (m1, n1)
        for (size_t m1_i = 0; m1_i < N; ++m1_i) {
            size_t n1_i = (m1_i + 1) % N;
            uint32_t m1 = tour[m1_i];
            uint32_t n1 = tour[n1_i];

            // for each (m2, n2) where m2 is choosen as the k-closest neighbor of m1
            for (size_t k = 0; k < nbhd.cols(); ++k) {
                uint32_t m2_i = whichSlot[nbhd.at(m1, k)];
                size_t n2_i = (m2_i + 1) % N;
                uint32_t m2 = tour[m2_i];
                uint32_t n2 = tour[n2_i];

                // if lower bound on new length is greater than upper bound of old length
                if (d.at(m1, m2) + minLink > d.at(m1, n1) + maxLink) {
                    break; // go to next m1
                }
                if (d.at(m1, m2) + d.at(n1, n2) < d.at(m1, n1) + d.at(m2, n2)) {
                    // make swap
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
inline void sortInTourOrder(vector<uint32_t>& cities, vector<size_t>& ind) {
    if (ind[0] < ind[4] && ind[0] > ind[2] 
        || ind[0] > ind[2] && ind[2] > ind[4]
        || ind[0] < ind[4] && ind[2] > ind[4]) {
            vector<uint32_t> tmp_c = cities;
            cities[0] = tmp_c[2]; cities[1] = tmp_c[3];
            cities[2] = tmp_c[0]; cities[3] = tmp_c[1];
            // cities[4] and cities[5] do not change
            vector<size_t> tmp_i = ind;
            ind[0] = tmp_i[2]; ind[1] = tmp_i[3];
            ind[2] = tmp_i[0]; ind[3] = tmp_i[1];
            // and their indices neither
        }
}

/**
 * Performs fast the 3-opt local optimization algorithm from the paper
 * "Large-Step Markov Chains for the Traveling Salesman Problem"
 */
inline void fastThreeOpt(Matrix& d, Matrix& nbhd, 
                     vector<uint32_t>& tour, vector<uint32_t>& whichSlot, 
                     uint32_t minLink, uint32_t& maxLink, 
                     chrono::steady_clock::time_point startTime, uint32_t timeLimit) {
    const size_t N = d.rows();

    bool improved = true;
    while (improved) {
        improved = false;

        if (timeOver(startTime, timeLimit))
            return; 

        // for each edge (m1, n1)
        for (size_t n1_i = 0; n1_i < N; ++n1_i) {
            size_t m1_i = (n1_i - 1 + N) % N;
            uint32_t m1 = tour[m1_i];
            uint32_t n1 = tour[n1_i];

            // for each edge (m2, n2) where n2 is choosen as the  k-nearest neighbor of m1
            for (size_t k1 = 0; k1 < nbhd.cols(); ++k1) {
                uint32_t n2_i = whichSlot[nbhd.at(m1, k1)];
                size_t m2_i = (n2_i + N - 1) % N; 
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
                for (size_t k2 = 0; k2 < nbhd.cols(); ++k2) {
                    uint32_t n3_i = whichSlot[nbhd.at(m1, k2)];
                    size_t m3_i = (n3_i + N - 1) % N;
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
                    vector<size_t> I = {m1_i, n1_i, m2_i, n2_i, m3_i, n3_i};
                    sortInTourOrder(C, I);

                    uint32_t oldDelta = d.at(C[0], C[1]) + d.at(C[2], C[3]) + d.at(C[4], C[5]);

                    // Try the two topologically different moves
                    // which correspond to four different ordered cases 
                    if (d.at(C[1], C[4]) + d.at(C[3], C[0]) + d.at(C[5], C[2]) < oldDelta) {
                        newBestMove = true;
                        reverseTourSegment(tour, whichSlot, I[0], I[5]);
                        reverseTourSegment(tour, whichSlot, I[1], I[2]);
                        reverseTourSegment(tour, whichSlot, I[3], I[4]);
                        maxLink = max(maxLink, max(max(d.at(C[1], C[4]), d.at(C[3], C[0])), d.at(C[5], C[2])));
                    } else if (d.at(C[3], C[0]) + d.at(C[5], C[1]) + d.at(C[2], C[4]) < oldDelta) {
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


int main(int argc, char *argv[]) {
    chrono::steady_clock::time_point startTime = chrono::steady_clock::now();
    const uint32_t timeLimit = 1900;
    const uint32_t twoOptTimeLimit = 200;
    const size_t K_NEAREST = 15;

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

    fastTwoOpt(distanceMatrix, nbhd, tour, whichSlot, minLink, maxLink, startTime, twoOptTimeLimit);
    fastThreeOpt(distanceMatrix, nbhd, tour, whichSlot, minLink, maxLink, startTime, timeLimit);

    for (int i = 0; i < tour.size(); ++i) {
        cout << tour[i] << endl;
    }
#if TESTING
    cout << "Tourlength: " << tourLength(tour, distanceMatrix) << endl;
#endif

    return 0;
}
