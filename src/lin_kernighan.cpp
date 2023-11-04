#include <set>
#include <cassert>
#include <iostream>
#include <list>
#include "lin_kernighan.h"

LK::LK(Tour& base, Matrix &distances, Neighbours &neighbours)
        : distances(distances), neighbours(neighbours), tour1(distances.dim()), tour2(distances.dim()), edges(), visited() {
    for (int i = 0; i < distances.dim(); ++i) {
        tour1[i] = base[i];
        tour2[i] = base[i];
    }
}

Tour &LK::get_tour() {
    if (iteration == 0) {
        return tour1;
    }
    return tour2;
}

bool LK::naive() {

    int n = distances.dim();

    std::vector<edge_t> X(n);
    std::vector<edge_t> Y(n);
    std::vector<length_t> Gi(n);
    std::vector<int> nodes(n);

    for (int i = 0; i < n; ++i) {
        X[i] = edge(0, 0);
        Y[i] = edge(-1, -1);
        Gi[i] = 0.0;
        nodes[i] = -1;
    }

    Tour& tour = iteration == 0 ? tour1 : tour2;
    Tour& T = iteration == 0 ? tour2 : tour1;

    int t1_idx = 0;
    if (!new_tour(X, Y, edge(-1, -1), 0)) {
        std::cerr << "Heuristic / Tour generation returns False" << std::endl;
    }

    Step2:
    int i = 1;
    // choose t1
    int t1 = tour[t1_idx];
    int x1_idx = 0;
    // TODO whenever we backtrack we have to reset to the last last

    Step3:
    // choose x1
    int t2 = tour.adjacent(t1_idx)[x1_idx];
    X[0] = edge(t1, t2);
    length_t gain = distances.at(t1, t2);

    std::list<std::pair<int, length_t>> potential_t3;

    for (auto t3 : neighbours.at(t2)) {
        if (tour.are_connected(t2, t3)) {
            continue;
        }
        length_t gi = gain - distances.at(t2, t3);
        if (gi > 0) {
            potential_t3.emplace_back(t3, gi);
        }
    }
    
    Step4:
    // choose y1 = (t2, t3)
    int t3 = -1;
    std::list<int> potential_y2;
    std::list<int> potential_x2;

    if (potential_t3.empty()) {
        goto Step12;
    }

    t3 = potential_t3.front().first;
    Gi[0] = potential_t3.front().second;
    potential_t3.pop_front();

    Y[0] = edge(t2, t3);
    // nodes[0] = t1
    // nodes[1] = t2
    nodes[2] = t3;

    // TODO populate potential_x2
    for (auto t4 : tour.adjacent(tour.index_of(t3))) {
        auto x2 = edge(t3, t4);
        // 6 (b)
        if (x2 == Y[0]) {
            continue;
        }
        // TODO 6 (a)
        potential_x2.push_back(t4);
    }

    Step5:
    i++;

    Step6:
    // xi = X[i - 1]

    if (i == 2) {
        if (potential_x2.empty()) {
            goto Step10;
        }
        auto even = potential_x2.front();
        X[1] = edge(nodes[2], even);
        potential_x2.pop_front();

        potential_y2.clear();

        for (auto y2 : neighbours.at(even)) {
            auto yi = edge(even, y2);
            if (tour.are_connected(even, y2)) {
                continue;
            }

            auto x1 = X[1];
            auto G2 = Gi[0] + distances.at(x1.first, x1.second) - distances.at(even, y2);
            // 7 (a)
            if (G2 <= 0.0) {
                continue;
            }

            // 7 (b)
            if (X[0] == yi || X[1] == yi) {
                continue;
            }

            // FIXME 7 (c)
            potential_y2.push_back(y2);
        }

        bool valid = new_tour(X, Y, edge(t1, even), 2);
        if (valid && T.length(distances) < tour.length(distances)) {
            std::cerr << T.length(distances) << " < " << tour.length(distances) << std::endl;
            iteration = (iteration + 1) % 2;
            return true;
        }
        nodes[3] = even;
    } else {
        // TODO choose first suitable xi that forms a tour
        bool found = false;
        int last = nodes[2 * i - 2];
        for (auto even : tour.adjacent(tour.index_of(last))) {
                X[i - 1] = edge(last, even);
            // 6 (b)
            bool fresh = true;
            for (int s = 0; s < i - 1; ++s) {
                if (Y[s] == X[i - 1]) {
                    fresh = false;
                    break;
                }
            }

            if (!fresh) {
                continue;
            }

            // 6 (a)
            bool valid = new_tour(X, Y, edge(t1, even), i);
            if (!valid) {
                continue;
            }

            if (T.length(distances) < tour.length(distances)) {
                std::cerr << T.length(distances) << " 2< " << tour.length(distances) << std::endl;
                iteration = (iteration + 1) % 2;
                return true;
            }

            found = true;
            nodes[2 * i - 1] = even;
            break;
        }
        if (!found) {
            goto Step8;
        }
        // TODO update last
    }

    Step7:
    if (i == 2) {
        if (potential_y2.empty()) {
            goto Step9;
        }

        auto odd = potential_y2.front();
        Y[1] = edge(nodes[3], odd);
        potential_y2.pop_front();
        auto x2 = X[1];
        Gi[1] = Gi[0] + distances.at(x2.first, x2.second) - distances.at(nodes[3], odd);

        nodes[4] = odd;
        goto Step5;
    } else {
        int last = nodes[2 * i - 1];
        for (auto odd : neighbours.at(last)) {
            auto yi = edge(last, odd);
            if (tour.are_connected(last, odd)) {
                continue;
            }
            auto xi = X[i - 1];
            Gi[i - 1] = Gi[i - 2] + distances.at(xi.first, xi.second) - distances.at(last, odd);
            // 7 (a)
            if (Gi[i - 1] <= 0.0) {
                continue;
            }

            // 7 (b)
            bool fresh = true;
            for (int s = 0; s < i; ++s) {
                if (X[s] == yi) {
                    fresh = false;
                    break;
                }
            }
            if (!fresh) {
                continue;
            }

            // FIXME 7 (c)

            Y[i - 1] = yi;
            nodes[2 * i] = odd;
            goto Step5;
        }
    }

    Step8:
    if (!potential_y2.empty()) {
        i = 2;
        goto Step7;
    }

    Step9:
    if (!potential_x2.empty()) {
        i = 2;
        goto Step6;
    }

    Step10:
    if (!potential_t3.empty()) {
        i = 1;
        goto Step4;
    }

    Step11:
    x1_idx++;
    if (x1_idx < 2) {
        i = 1;
        goto Step3;
    }
    
    Step12:
    t1_idx++;
    if (t1_idx < n) {
        goto Step2;
    }

    std::cerr << "okay" << std::endl;
    std::cerr << T.length(distances) << " " << T.is_valid() << std::endl;
    std::cerr << tour.length(distances) << std::endl;
    return false;
}

bool LK::new_tour(std::vector<edge_t> &X, std::vector<edge_t> &Y, edge_t final, int i) {
    int n = distances.dim();

    Tour& tour = iteration == 0 ? tour1 : tour2;
    Tour& T = iteration == 0 ? tour2 : tour1;

    // used for validating the tour
    visited.clear();
    for (int j = 0; j < n; ++j) {
        visited.insert(j);
    }

    // keeps track of all edges in new tour
    // edges = (tour-edges UNION added) SETDIFF removed
    edges.clear();
    edges.insert(final);
    for (int j = 0; j < n; ++j) {
        edge_t e = edge(tour[j], tour[(j + 1) % n]);
        bool remove = false;
        for (int k = 0; k < i; ++k) {
            if (X[k] == e) {
                remove = true;
                break;
            }
        }
        if (!remove) {
            edges.insert(e);
        }
    }
    for (int j = 0; j < i - 1; ++j) {
        if (Y[j] == edge(-1, -1) || Y[j].first == Y[j].second) {
            std::cerr << "ERR" << std::endl;
        }
        edges.insert(Y[j]);
    }

    // FIXME maybe this does not have to hold
    //if (edges.size() != n) {
      //  return false;
    //}

    int idx = 0;
    // always start at 0
    T[0] = idx;

    for (int j = 1; j < n; ++j) {
        edge_t found_edge(-1, -1);
        bool not_found = true;
        // search for some edge that has idx as an endpoint and choose it
        for (auto& e : edges) {
            if (e.first == idx) {
                idx = e.second;
                found_edge = e;
                not_found = false;
                break;
            } else if (e.second == idx) {
                idx = e.first;
                found_edge = e;
                not_found = false;
                break;
            }
        }
        if (not_found) {
            return false;
        }
        T[j] = idx;
        // remove found edge to (a) increase performance, and (b) don't loop e.g. 1 -> 2, 2 <- 1, 1 -> 2, 2 <- 1, ...
        if (!not_found) {
            edges.erase(found_edge);
        }
    }

    for (int j = 0; j < n; ++j) {
        visited.erase(T[j]);
    }

    return visited.empty();
}
