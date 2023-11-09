#include "fastOptTSP.cpp"
#include "Minimum-Cost-Perfect-Matching/Graph.h"
// steps:
// 1. construct a minimal spanning tree (linear?)
// 2. for odd degree nodes in that tree: find minimal weight (distance) perfect matching (O(n^3))
//     note: assignment pages says we are allowed to use an existing implementation if we cite the source
//     for example: https://github.com/dilsonpereira/Minimum-Cost-Perfect-Matching
// 3. combine matching and tree into one graph, now all nodes are even degree (linear?)
// 4. construct Euler tour of that new graph (linear time))
// 5. TSP tour is the Euler tour with all duplicates after first occurence removed
//
// uses https://github.com/dilsonpereira/Minimum-Cost-Perfect-Matching
// and its accompanying graph structure (undirected graphs)

using namespace std;

pair< Graph, vector<uint32_t> > matrix_to_graph(Matrix& d) {
  int n = d.rows();
  Graph g(n);
  vector<uint32_t> cost(n*(n-1)/2);

  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      g.AddEdge(i, j);
      cost[g.GetEdgeIndex(i, j)] = d.at(i, j);
    }
  }
  return make_pair(g, cost);
}

// create min weight spanning tree (Prim's algorithm)
// can be adapted into a linear time algo that uses Prim's
// https://en.wikipedia.org/wiki/Prim%27s_algorithm
pair<Graph, vector<uint32_t>> min_spanning_tree(Graph& g, vector<uint32_t>& cost) {
  int n = g.GetNumVertices();
  vector<uint32_t> cheapest(n, numeric_limits<uint32_t>::max());
  // edge indices in g
  vector<uint32_t> cheapest_idx(n, -1);
  Graph forest(n);
  bool not_added[n];
  for (int i = 0; i < n; i++) {
    not_added[i] = true;
  }
  int num_not_added = n;
  //process ("add to forest") one vertex at a time to build the full tree
  for (int b = 0; b < n; b++) {
    uint32_t min_val = numeric_limits<uint32_t>::max();
    int min_vert = 0;
    //find unadded vertex with (as far as we know) cheapest edge
    for (int i = 0; i < n; i++) {
      if (cheapest[not_added[i]] < min_val) {
        min_val = cheapest[not_added[i]];
        min_vert = i;
      }
    }
    //count as added
    not_added[min_vert] = false;
    for (int i = 0; i < n; i++) {
      if (i != min_vert) {
        if (not_added[i] && cost[g.GetEdgeIndex(min_vert, i)] < cheapest[i]) {
          cheapest[i] = cost[g.GetEdgeIndex(min_vert, i)];
          cheapest_idx[i] = g.GetEdgeIndex(min_vert, i);
        }
      }
    }
  }
  // add edges to forest
  vector<uint32_t> cost_forest(n-1);
  for (int i = 0; i < n; i++) {
    // should be all vertices except the first
    if (cheapest_idx[i] != -1) {
      pair<int, int> edge = g.GetEdge(cheapest_idx[i]);
      forest.AddEdge(edge.first, edge.second);
      cost_forest[forest.GetEdgeIndex(edge.first, edge.second)] = cost[cheapest_idx[i]];
    }
  }
  return make_pair(forest, cost_forest);
}

// subgraph with only the odd degree nodes
// NOT READY is this supposed to be all edges from the graph, or just those in the tree??
// started to implement for taking only the edges that are in the spanning tree, but idk
pair< Graph, vector<uint32_t> > odd_nodes_graph(Graph& g, vector<uint32_t>& cost) {
  int n = g.GetNumVertices();

  int odd_map[n];
  int counter = 0;
  for(int i = 0; i < n; i++) {
    if (g.AdjList(i).size() % 2 == 1) {
      odd_map[i] = counter;
      counter++;
    } else {
      odd_map[i] = -1;
    }
  }

  Graph o(counter);
  vector<uint32_t> cost_odd;
  // now the vertices don't have the same numbers anymore, what do?
  // return a map or find a way to keep all verts in the graph?

  for(int i = 0; i < n; i++) {
    if (odd_map[i] != -1) {
      list<int> adj_list = g.AdjList(i);
      for (std::list<int>::iterator it = adj_list.begin(); it != adj_list.end(); it++) {
        int e = odd_map[*it];
        // check if edge goes to another odd vert, and if so add edge
        if (e != -1) {
          o.AddEdge(odd_map[i], e);
          cost_odd[o.GetEdgeIndex(odd_map[i], e)] = cost[g.GetEdgeIndex(i, *it)];
        }
      }
    }
  }

  
  
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

  vector<uint32_t> tour(n);
  pair<Graph, vector<uint32_t>> converted = matrix_to_graph(d);
  Graph g = converted.first;
  vector<uint32_t> cost = converted.second;
  pair<Graph, vector<uint32_t>> min_span_tree = min_spanning_tree(g, cost);
  pair<Graph, vector<uint32_t>> odd = odd_nodes_graph(min_span_tree.first, min_span_tree.second);
  odd_nodes(tree, odd);
  Matrix matching(n, n);
  min_perfect_matching(odd);
  Matrix out(n, n);
  combine_graphs(odd, tree, out);
  euler_tour(out, tour);
  remove_duplicates(tour);
  
  return tour;
  
}
