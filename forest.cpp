#include <RcppArmadillo.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;

using namespace Rcpp;
using namespace arma;

// declaration for new node
class Node {
 public:
  int idx;
  vector<int> adjList;

  Node(){};

  Node(int x) { idx = x; }

  // int getIdx() { return idx; };
  // void setIdx(int idx_) { idx = idx_; };

  bool checkIfConnectTO(int b) {
    return std::find(adjList.begin(), adjList.end(), b) != adjList.end();
  }

  void connectTo(int b, map<int, Node> &graph) {
    // check if adjList includes b already
    if (!checkIfConnectTO(b)) {
      adjList.push_back(b);
      // add this node's idx to the graph[b] too
      if (!graph[b].checkIfConnectTO(idx)) graph[b].adjList.push_back(idx);
    }
  }

  void cutFrom(int b, map<int, Node> &graph) {
    if (checkIfConnectTO(b)) {
      adjList.erase(std::remove(adjList.begin(), adjList.end(), b),
                    adjList.end());
      graph[b].cutFrom(idx, graph);
    }
  }

  void getConnectedNodesAux(int parent_idx, vector<int> &node_list,
                            map<int, Node> graph, int max_num_neighbors_to_check) {
      
    if(node_list.size()<max_num_neighbors_to_check){
        node_list.emplace_back(idx);

        for (int i = 0; i < adjList.size(); ++i) {
          if (adjList[i] != parent_idx) {
            graph[adjList[i]].getConnectedNodesAux(idx, node_list, graph, max_num_neighbors_to_check);
          }
        }
    }
  }

  vector<int> getConnectedNodes(map<int, Node> graph, int max_num_neighbors_to_check) {
    vector<int> node_list;

//     assert(max_num_neighbors_to_check >0);
      
    getConnectedNodesAux(idx, node_list, graph,max_num_neighbors_to_check);

    return node_list;
  }
};

class Graph {
 public:
  map<int, Node> graph;
  int n;
  mat S;
  mat A;

  vector<pair<int, int>> edge_list;

  Graph(int n_) {
    n = n_;
    for (int i = 0; i < n; i++) {
      graph[i] = Node(i);
    }
    A.zeros(n, n);
  };

  void connect(int a, int b, bool push_to_list = false) {
    if (a < n and b < n and a != b) {
      graph[a].connectTo(b, graph);
      // A(a, b) = 1;
      // A(b, a) = 1;

      if (push_to_list) {
        edge_list.push_back(make_pair(a, b));
      }
    }
  };

  void cut(int a, int b) {
    if (a < n and b < n) {
      graph[a].cutFrom(b, graph);
      // A(a, b) = 0;
      // A(b, a) = 0;
    }
  }

  IntegerVector getConnectedNodes(int idx) {
    IntegerVector r_node_list;

    if (idx < n) {
      auto node_list = graph[idx].getConnectedNodes(graph,n);
      for (int i : node_list) r_node_list.push_back(i);
    }
    return r_node_list;
  }

  arma::uvec arma_setdiff(arma::uvec x, arma::uvec y) {
    std::vector<int> a = arma::conv_to<std::vector<int>>::from(arma::sort(x));
    std::vector<int> b = arma::conv_to<std::vector<int>>::from(arma::sort(y));
    std::vector<int> out;

    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(out, out.end()));

    return arma::conv_to<arma::uvec>::from(out);
  }

  uvec getConnectedNodesArma(int idx, int max_num_neighbors_to_check) {
    uvec arma_node_list;

    if (idx < n) {
      auto node_list = graph[idx].getConnectedNodes(graph, max_num_neighbors_to_check);
      arma_node_list = conv_to<uvec>::from(node_list);
    }
    return arma_node_list;
  }

  void setS(NumericMatrix S_) { S = as<arma::mat>(S_); }

  NumericMatrix getA() {
    A.zeros(n, n);
    for (int k = 0; k < edge_list.size(); k++) {
      int a = edge_list[k].first;
      int b = edge_list[k].second;
      A(a, b) = 1;
      A(b, a) = 1;
    }
    return wrap(A);
  }

  void updateGraph(IntegerVector &C_r, double lam, int max_num_neighbors_to_check) {
    uvec C = conv_to<uvec>::from(as<std::vector<int>>(C_r));

    uvec full_idx = linspace<uvec>(0, n - 1, n);
    uvec one_n_1 = linspace<uvec>(1, n - 1, n - 1);

    for (int k = 0; k < edge_list.size(); k++) {
      pair<int, int> edge = edge_list[k];
      //          Rcout<< edge.first<<endl;
      //          Rcout<< edge.second<<endl;

      cut(edge.first, edge.second);

      uvec Ga = getConnectedNodesArma(edge.first, max_num_neighbors_to_check);
      uvec Gb = getConnectedNodesArma(edge.second, max_num_neighbors_to_check);

      if (sum(Gb == 0) >0) {
        // swap Ga and Gb
        uvec temp = Ga;
        Ga = Gb;
        Gb = temp;
      }
        
//       Rcout<< Ga.n_elem<<" "<< Gb.n_elem<<endl;

      // assert(sum(Ga == 0) > 0);
      // assert(sum(Gb == 0) == 0);

      int Kstar = 1;

      uvec uniqueCinGa = unique(C(Ga));
      int maxCGa = max(uniqueCinGa);

      if (Ga.n_elem > 1) {
        Kstar = uniqueCinGa.n_elem - 1;
      }
        

      if (Ga.n_elem > 0 && Gb.n_elem > 0) {
        mat S_GaGb = S.submat(Ga, Gb);

        mat gumbel = -log(-log(randu(Ga.n_elem, Gb.n_elem)));
        S_GaGb = S_GaGb + gumbel;
          
        if (sum(Ga == 0) >0) { // if 0 is in Ga
            S_GaGb.row(0) = S_GaGb.row(0) + log(lam) - log((double)Kstar + 1);
        }
          
        uword max_idx = S_GaGb.index_max();

        uvec idx2 = ind2sub(size(S_GaGb), max_idx);

        int new_i =  Ga(idx2(0));
        int new_j =  Gb(idx2(1));

        connect(new_i, new_j, false);
        edge_list[k] = make_pair(new_i, new_j);

        if (new_i == 0) {  // add a cluster
          uvec one_n_1 = linspace<uvec>(1, maxCGa+1, maxCGa+1);
          int new_C_idx =  min(arma_setdiff(one_n_1, uniqueCinGa));
          //             Rcout<<new_C_idx<<endl;
          C(Gb).fill(new_C_idx);
        } else {
          C(Gb).fill(C(new_i));
        }
      } else {
        // restore the link
        connect(edge.first, edge.second, false);
      }
    
        
        
    }
      
    for (int i = 0; i < C.n_elem; i++) {
      C_r[i] = C(i);
    }
  }
};

RCPP_MODULE(mod) {
  class_<Graph>("Graph")

      .constructor<int>()

      .method("connect", &Graph::connect)
      .method("cut", &Graph::cut)
      .method("getA", &Graph::getA)
      .method("setS", &Graph::setS)
      .method("updateGraph", &Graph::updateGraph)

      .method("getConnectedNodes", &Graph::getConnectedNodes);
}
