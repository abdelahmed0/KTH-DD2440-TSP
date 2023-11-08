#include "src/fastOptTSP.cpp"
// steps:
// 1. construct a minimal spanning tree (linear?)
// 2. for odd degree nodes in that tree: find minimal weight (distance) perfect matching (O(n^3))
//     note: assignment pages says we are allowed to use an existing implementation if we cite the source
//     for example: https://github.com/dilsonpereira/Minimum-Cost-Perfect-Matching
// 3. combine matching and tree into one graph, now all nodes are even degree (linear?)
// 4. construct Euler tour of that new graph (linear time))
// 5. TSP tour is the Euler tour with all duplicates after first occurence removed

using namespace std;

// create min weigh spanning tree
void min_spanning_tree(Matrix& d, Matrix& tree) {
  
}

// new matrix where all edges not between odd nodes = inf
void odd_nodes(Matrix& m, Matrix& odd) {
  
}

// mutate in place
void min_perfect_matching(Matrix& d) {
  
}

void combine_graphs(Matrix& a, Matrix& b, Matrix& out) {
  
}

void euler_tour(Matrix& m, vector<uint32_t>& tour) {
}

vector<uint32_t> remove_duplicates(vector<uint32_t>& tour) {
  return tour;
}

vector<uint32_t> christofides(Matrix& d) {

  int n = d.rows();
  vector<uint32_t> tour(n);
  Matrix tree(n, n);
  min_spanning_tree(d, tree);
  Matrix odd(n, n);
  odd_nodes(tree, odd);
  Matrix matching(n, n);
  min_perfect_matching(odd);
  Matrix out(n, n);
  combine_graphs(odd, tree, out);
  euler_tour(out, tour);
  remove_duplicates(tour);
  
  return tour;
  
}
