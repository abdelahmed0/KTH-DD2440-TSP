#include <set>
#include <cassert>
#include <iostream>
#include "lin_kernighan.h"

LK::LK(Tour &tour, Matrix &distances, Neighbours &neighbours)
        : distances(distances), neighbours(neighbours), tour(tour), out(distances.dim()) {
    for (int i = 0; i < distances.dim(); ++i) {
        out[i] = tour[i];
    }
}

Tour &LK::get_tour() {
    return tour;
}

bool LK::step() {
    // idxi correspond to indices in the tour array
    // ti correspond to city index i.e. element at idxi
    // ti = tour[idxi]

    // we implicitly keep track of the current iteration by keeping the edge set, i.e., i = removed_set.size()
    int n = distances.dim();
    auto curr_length = tour.length(distances);
    for (int idx1 = 0; idx1 < n; ++idx1) { // Step 2
        auto t1 = tour[idx1];
        for (auto t2 : {tour.successor(idx1)}) { // Step 3
            // considering adjacent is equivalent to always picking edges \in T
            auto x1 = edge(t1, t2);
            auto cx1 = distances.at(t1, t2);

            // FIXME: maybe both X and Y can be moved into the 6/7 loop
            std::set<edge_t> X1; // X = removed
            X1.insert(x1);

            for (auto t3 : neighbours.at(t2)) { // Step 4
                auto idx3 = tour.index_of(t3);

                // this checks (t2, t3) \notin T
                if (tour.are_connected(t2, t3)) {
                    continue;
                }

                auto y1 = edge(t2, t3);
                auto cy1 = distances.at(t2, t3);

                auto gain = cx1 - cy1;
                std::set<edge_t> Y1; // Y = added
                Y1.insert(y1);

                if (gain > 0) {
                    // Step 6
                    // choose x_i = (t_2i-1, t_2i) \in T such that
                    // (a) if t_2i is joined to t_1, the resulting config is a tour, T', and
                    // (b) x_i != y_s for all s < i (NOTE: this can be done by checking that x_i \notin Y as we have an iterative process)
                    // If T' is a better tour than T, let T = T' and go to Step 2 (NOTE: in our case we could exit this loop to have some more control)

                    // we split step 6/7 into two parts either i = 2 or not i > 2

                    // Step 6.1: i = 2
                    for (auto t4 : tour.adjacent(idx3)) {
                        auto X2(X1);

                        auto x2 = edge(t3, t4);
                        auto cx2 = distances.at(t3, t4);
                        // (b)
                        if (Y1.count(x2) > 0) {
                            continue;
                        }

                        X2.insert(x2);
                        auto Y2_tmp(Y1);
                        Y2_tmp.insert(edge(t4, t1));
                        Tour T(n);
                        // we ignore (a) for i = 2
                        bool valid = generate_new_tour(T, Y2_tmp, X2);
                        if (valid && T.length(distances) < curr_length) {
                            tour = T;
                            return true;
                        }

                        // Step 7
                        // choose y_i = (t_2i, t_2i+1) \notin T such that
                        // (a) G_i > 0,
                        // (b) y_i != x_s for all s <= i, and
                        // (c) x_i+1 exists
                        // If such y_i exists, go to Step 5

                        // Step 7.1: i = 2
                        for (auto t5 : neighbours.at(t4)) {
                            auto G2 = gain + (cx2 - distances.at(t4, t5));
                            auto y2 = edge(t4, t5);
                            // (a) & (b)
                            if (G2 <= 0.0 || X2.count(y2) > 0) {
                                continue;
                            }

                            auto Y2(Y1);
                            Y2.insert(y2);

                            // FIXME ignore condition (c) for now

                            bool success = recurse(T, X2, Y2, t1, t5, G2);
                            if (success) {
                                return true;
                            }
                        }
                    }

                    // Step 8
                    // if there is an untried alternative for y_2, let i = 2 and go to Step 7

                    // Step 9
                    // if there is an untried alternative for x_2, let i = 2 and go to Step 6
                }
                // repeat choice for step 4 and then try other t2
                // Step 10
            }
            // if both choices for t2 don't work, re-decide t1
            // Step 11
        }
        // try every node as t1 otherwise forfeit
        // Step 12
    }
    return false;
}

bool LK::recurse(Tour &T, std::set<edge_t> Xi, std::set<edge_t> Yi, int t1, int last, length_t Gi) {
    // Step 6
    for (auto t2i : tour.adjacent(tour.index_of(last))) {
        auto xi = edge(last, t2i);
        if (Yi.count(xi) > 0 || Xi.count(xi) > 0) {
            continue;
        }

        auto Y_tmp(Yi);
        auto X_tmp(Xi);
        Y_tmp.insert(edge(t2i, t1));
        X_tmp.insert(xi);

        bool valid = generate_new_tour(T, Y_tmp, X_tmp);
        if (!valid) {
            continue;
        }
        Xi.insert(xi);
        if (T.length(distances) < tour.length(distances)) {
            tour = T;
            return true;
        }

        // Step 7
        for (auto t2i1 : neighbours.at(t2i)) {
            auto yi = edge(t2i, t2i1);
            if (tour.are_connected(t2i, t2i1) || Xi.count(yi) > 0 || Yi.count(yi) > 0) {
                continue;
            }
            auto gi = distances.at(xi.first, xi.second) - distances.at(yi.first, yi.second);
            if (Gi + gi > 0) {
                auto tmp(Yi);
                tmp.insert(yi);
                // we have to cha
                if (recurse(T, Xi, tmp, t1, t2i1, Gi + gi)) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool LK::generate_new_tour(Tour &t2, std::set<edge_t> &added, std::set<edge_t> &removed) {
    int n = distances.dim();
    // used for validating the tour
    std::set<int> visited;

    // keeps track of all edges in new tour
    // edges = (tour-edges UNION added) SETDIFF removed
    std::set<edge_t> edges;
    for (int i = 0; i < n; ++i) {
        visited.insert(i);

        edge_t e = edge(tour[i], tour[(i + 1) % n]);
        if (!set_contains(removed, e)) {
            edges.insert(e);
        }
    }
    for (auto e : added) {
        edges.insert(e);
    }

    // FIXME maybe this does not have to hold
    if (edges.size() != n) {
        /*std::cerr << edges.size() << " |Y|=" << added.size() << " |X|=" << removed.size() << "\t Y= ";
        for (auto [e1, e2] : added) {
            std::cerr << "(" << e1 << ", "  << e2 << ") ";
        }
        std::cerr << "\tX= ";
        for (auto [e1, e2] : removed) {
            std::cerr << "(" << e1 << ", "  << e2 << ") ";
        }
        std::cerr << std::endl;*/
        return false;
    }

    int idx = 0;
    // always start at 0
    t2[0] = idx;

    // TODO check if works
    for (int i = 1; i < n; ++i) {
        edge_t found_edge(-1, -1);
        // search for some edge that has idx as an endpoint and choose it
        for (auto& e : edges) {
            if (e.first == idx) {
                idx = e.second;
                found_edge = e;
            } else if (e.second == idx) {
                idx = e.first;
                found_edge = e;
            }
        }
        t2[i] = idx;
        // remove found edge to (a) increase performance, and (b) don't loop e.g. 1 -> 2, 2 <- 1, 1 -> 2, 2 <- 1, ...
        if (found_edge != std::pair(-1, -1)) {
            edges.erase(found_edge);
        }
    }

    for (int i = 0; i < n; ++i) {
        visited.erase(t2[i]);
    }

    return visited.size() == n;
}