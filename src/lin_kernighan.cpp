#include <set>
#include "lin_kernighan.h"

LK::LK(Tour &tour, Matrix &distances, Neighbours &neighbours)
        : distances(distances), neighbours(neighbours), tour(tour), out(tour) {}


bool LK::step() {
    // idxi correspond to indices in the tour array
    // ti correspond to city index i.e. element at idxi
    // ti = tour[idxi]

    // we implicitly keep track of the current iteration by keeping the edge set
    int n = distances.dim();
    auto curr_length = tour.length(distances);
    for (int idx1 = 0; idx1 < n; ++idx1) { // Step 2
        auto t1 = tour[idx1];

        for (auto t2 : tour.adjacent(idx1)) { // Step 3
            // considering adjacent is equivalent to always picking edges \in T
            auto x1 = edge(t1, t2);
            auto cx1 = distances.at(t1, t2);
            // FIXME: maybe both X and Y can be moved into the 6/7 loop
            std::set<edge_t> removed_set; // X
            removed_set.insert(x1);

            for (auto t3 : neighbours.at(t2)) { // Step 4
                auto idx3 = tour.index_of(t3);
                auto adj = tour.adjacent(idx3);
                // checking that t2 is adjacent to t3
                // this checks (t2, t3) \notin T
                if (adj[0] == t2 || adj[1] == t2) {
                    continue;
                }

                auto y1 = edge(t2, t3);
                auto cy1 = distances.at(t2, t3);

                auto gain = cx1 - cy1;
                std::set<edge_t> added_set; // Y
                added_set.insert(y1);

                if (gain > 0) {
                    // Step 6
                    // choose x_i = (t_2i-1, t_2i) \in T such that
                    // (a) if t_2i is joined to t_1, the resulting config is a tour, T', and
                    // (b) x_i != y_s for all s < i (NOTE: this can be done by checking that x_i \notin Y as we have an iterative process)
                    // If T' is a better tour than T, let T = T' and go to Step 2 (NOTE: in our case we could exit this loop to have some more control)

                    auto t2i_1 = t3;
                    auto X(removed_set);
                    auto Y(added_set);

                Step_6:
                    auto idx2i_1 = tour.index_of(t3);
                    for (auto t2i : tour.adjacent(idx2i_1)) {
                        auto xi = edge(t2i_1, t2i);
                        // 6 (b)
                        if (added_set.count(xi) > 0) {
                            continue;
                        }

                        auto relink = edge(t2i, t1);

                        X.insert(xi);
                        Y.insert(relink);

                        Tour t_prime = tour.update(Y, X);
                        bool valid_tour = t_prime.is_valid();
                        // 6 (a)
                        if (!valid_tour && X.size() > 2) {
                            continue;
                        }

                        if (valid_tour && t_prime.length(distances) < curr_length) {
                            tour = t_prime;
                            return true;
                        }

                        for (auto t2i1 : neighbours.at(t2i)) {

                        }
                    }


                    // Step 7
                    // choose y_i = (t_2i, t_2i+1) \notin T such that
                    // (a) G_i > 0,
                    // (b) y_i != x_s for all s <= i, and
                    // (c) x_i+1 exists
                    // If such y_i exists, go to Step 5

                    // Step 8
                    // if there is an untried alternative for y_2, let i = 2 and go to Step 7

                    // Step 9
                    // if there is an untried alternative for x_2, let i = 2 and go to Step 6
                    return true;
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

bool LK::chooseX(int t1, int t2, length_t& gain, std::set<edge_t>& removed_edges, std::set<edge_t>& added_edges) {
    for (auto t2i : tour.adjacent(t2)) {
        auto xi = edge(t2, t2i);
        auto gi = gain + distances.at(t2, t2i);
        if (!set_contains(added_edges, xi) && !set_contains(removed_edges, xi)) {
            std::set<edge_t> add(added_edges);
            std::set<edge_t> rem(removed_edges);
            rem.insert(xi);

            auto xr = edge(t2i, t1);
            add.insert(xr);

            auto relink = gi - distances.at(t2i, t1);
            Tour new_tour = tour.update(add, rem);
            if (!new_tour.is_valid() && add.size() > 2) { // FIXME
                continue;
            }

            if (new_tour.is_valid() && relink > 0) {
                out = new_tour;
                return true;
            } else {
                bool choice = chooseY(t1, t2i, gi, rem, added_edges);
                if (choice && removed_edges.size() == 2) {
                    return true;
                }
                return choice;
            }
        }
    }
    return false;
}

bool LK::chooseY(int t1, int t2, length_t& gain, std::set<edge_t>& removed_edges, std::set<edge_t>& added_edges) {
    for (auto node : neighbours.at(t2)) {
        auto yi = edge(t2, node);
        auto add(added_edges);
        add.insert(yi);

        auto gi = gain - distances.at(t2, node);
        if (chooseX(t1, node, gi, removed_edges, add)) {
            return true;
        }

    }
    return false;
}

Tour &LK::get_tour() {
    return out;
}
