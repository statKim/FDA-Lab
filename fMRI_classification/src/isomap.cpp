#include <Rcpp.h>
//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <limits>
#include "trapzRcpp.h"

using namespace Rcpp;
using namespace Eigen;
using namespace std;



// L2 distance between x and y from timepoints t
// [[Rcpp::export]]
double dist_L2(NumericVector x, NumericVector y, NumericVector t) {
  return sqrt(trapzRcpp(t, pow(x - y, 2)));
}


// Distance matrix using L2 distance
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_mat_L2(Rcpp::NumericMatrix X, NumericVector t) {
  int n = X.nrow();
  // int p = X.ncol();
  
  // Rcpp::Rcout << n << " " << p << "\n";
  
  Rcpp::NumericMatrix dist_mat(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i > j) {
        dist_mat(i, j) = dist_L2(X.row(i), X.row(j), t);
        dist_mat(j, i) = dist_mat(i, j);
      }
    }
  }
  // n <- nrow(X)
  // p <- ncol(X)
  // 
  // dist_mat <- matrix(0, n, n)
  // for (i in 1:n) {
  //   for (j in 1:n) {
  //     if (i > j) {
  //       dist_mat[i, j] <- dist_l2(X[i, ], X[j, ], t)
  //     }
  //   }
  // }
  
  return dist_mat;
}



/// Dijkstra algorithm from ChatGPT 3.5
#include <Rcpp.h>
#include <vector>
#include <limits>

using namespace Rcpp;
using namespace std;

// Function to find the vertex with minimum distance value, from the set of vertices
// not yet included in shortest path tree
int minDistance(vector<double>& dist, vector<bool>& sptSet, int V) {
  // Initialize min value
  double min = numeric_limits<double>::max();
  int min_index = -1; // Initialize min_index to avoid uninitialized value
  
  for (int v = 0; v < V; v++) {
    if (sptSet[v] == false && dist[v] <= min) {
      min = dist[v];
      min_index = v;
    }
  }
  return min_index;
}

// Dijkstra algorithm function
// [[Rcpp::export]]
List dijkstra_cpp(NumericMatrix graph, int src) {
  int V = graph.nrow();
  vector<double> dist(V);  // The output array. dist[i] will hold the shortest distance from src to i
  vector<int> parent(V);   // Store the parent node of each node in the shortest path
  
  // sptSet[i] will be true if vertex i is included in shortest
  // path tree or shortest distance from src to i is finalized
  vector<bool> sptSet(V, false);
  
  // Initialize all distances as INFINITE and src as source
  for (int i = 0; i < V; i++) {
    dist[i] = numeric_limits<double>::max();
    parent[i] = -1; // Initialize parent array
  }
  dist[src] = 0;  // Distance of source vertex from itself is always 0
  
  // Find shortest path for all vertices
  for (int count = 0; count < V - 1; count++) {
    // Pick the minimum distance vertex from the set of vertices
    // not yet processed. u is always equal to src in the first
    // iteration.
    int u = minDistance(dist, sptSet, V);
    
    // Mark the picked vertex as processed
    sptSet[u] = true;
    
    // Update dist value of the adjacent vertices of the
    // picked vertex.
    for (int v = 0; v < V; v++) {
      // Update dist[v] only if is not in sptSet, there is an
      // edge from u to v, and total weight of path from src to
      // v through u is smaller than current value of dist[v]
      if (!sptSet[v] && graph(u,v) && dist[u] != numeric_limits<double>::max() && dist[u] + graph(u,v) < dist[v]) {
        dist[v] = dist[u] + graph(u,v);
        parent[v] = u; // Update parent node
      }
    }
  }
  
  // Construct the paths
  List paths(V);
  for (int i = 0; i < V; i++) {
    int node = i;
    IntegerVector path;
    while (node != src && parent[node] != -1) { // Modify this line
      path.push_front(node);
      node = parent[node];
    }
    path.push_front(src);
    paths[i] = path;
  }
  
  NumericVector dist_output(V);
  for(int i = 0; i < V; i++)
    dist_output[i] = dist[i];
  
  return List::create(Named("distances") = dist_output,
                      Named("paths") = paths);
}


// Compute actual path distance after Dijkstra algorithm
// [[Rcpp::export]]
NumericVector find_path_dist(NumericMatrix dist_mat, List path) {
  int n = dist_mat.nrow();
  int m;
  
  NumericVector path_dist(n);  // zero vector
  for (int i = 0; i < n; i++) {
    NumericVector path_i = path[i];
    m = path_i.size();
    for (int j = 1; j < m; j++) {
      path_dist[i] = path_dist[i] + dist_mat(path_i[j-1], path_i[j]);
    }
  }
  
  return path_dist;
}

// Compute actual path distance matrix after Dijkstra algorithm
// - dist_mat_p: penalized distance matrix
// - dist_mat: original distance matrix
// [[Rcpp::export]]
NumericMatrix geod_dist_mat(NumericMatrix dist_mat_p, NumericMatrix dist_mat) {
  int n = dist_mat.nrow();
  List dijkstra_res;
  
  NumericMatrix geo_path_mat(n, n);  // zero matrix
  for (int i = 0; i < n; i++) {
    dijkstra_res = dijkstra_cpp(dist_mat_p, i);
    geo_path_mat.row(i) = find_path_dist(dist_mat, dijkstra_res["paths"]);
  }
  
  return geo_path_mat;
}


// 
// /*** R
// d <- matrix(c(0,Inf,Inf,Inf,Inf,7,0,2,3,Inf,4,Inf,0,Inf,Inf,6,Inf,5,0,1,1,Inf,Inf,Inf,0), 5, 5)
// d <- matrix(c(0,0,0,0,0,7,0,2,3,0,4,0,0,0,0,6,0,5,0,1,1,0,0,0,0), 5, 5)
// d <- d + t(d)
// # dijkstra(d, 1)
// dijkstra_res <- dijkstra_cpp(d, 0)
// # 0 5 4 2 1
// geod_dist_mat(d, d)
// */

// 
// /*** R
// # gr <- seq(0, 1, length.out = 10)
// # X <- matrix(rnorm(200), 20, 10)
// # as.dist(est_geo_dist(X, gr, eps=0, delta=0))
// # 예제 그래프
// graph <- matrix(c(
//   0, 4, 0, 0, 0, 0, 0, 8, 0,
//   4, 0, 8, 0, 0, 0, 0, 11, 0,
//   0, 8, 0, 7, 0, 4, 0, 0, 2,
//   0, 0, 7, 0, 9, 14, 0, 0, 0,
//   0, 0, 0, 9, 0, 10, 0, 0, 0,
//   0, 0, 4, 14, 10, 0, 2, 0, 0,
//   0, 0, 0, 0, 0, 2, 0, 1, 6,
//   8, 11, 0, 0, 0, 0, 1, 0, 7,
//   0, 0, 2, 0, 0, 0, 6, 7, 0), 
//   nrow=9, byrow=TRUE)
// 
// eps <- 5
// graph <- ifelse(graph > eps, Inf, graph)
// 
// # Dijkstra 알고리즘 호출
// dist_mat <- matrix(0, nrow(graph), ncol(graph))
// for (i in 1:nrow(graph)) {
//   dist_mat[i, ] <- dijkstra_cpp(graph, i-1)
// }
// dist_mat
// */
