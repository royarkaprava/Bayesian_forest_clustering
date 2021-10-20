#include <RcppArmadillo.h>
#include <sstream>

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;

using namespace Rcpp;

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
      if (!graph[b].checkIfConnectTO(idx))
        graph[b].adjList.push_back(idx);
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
                            map<int, Node> graph) {

    node_list.push_back(idx);

    for (int i = 0; i < adjList.size(); ++i) {
      if (adjList[i] != parent_idx) {
        graph[adjList[i]].getConnectedNodesAux(idx, node_list, graph);
      }
    }
  }

  vector<int> getConnectedNodes(map<int, Node> graph) {
    vector<int> node_list;

    getConnectedNodesAux(idx, node_list, graph);

    return node_list;
  }
};

class Graph {
public:
  map<int, Node> graph;
  int n;

  Graph(int n_) {
    n = n_;
    for (int i = 0; i < n; i++) {
      graph[i] = Node(i);
    }
  };

  void connect(int a, int b) {
    if (a < n and b < n) {
      graph[a].connectTo(b, graph);
    }
  };

  void cut(int a, int b) {
    if (a < n and b < n) {
      graph[a].cutFrom(b, graph);
    }
  }

  IntegerVector getConnectedNodes(int idx) {

    IntegerVector r_node_list;

    if (idx < n) {
      auto node_list = graph[idx].getConnectedNodes(graph);
      for (int i : node_list)
        r_node_list.push_back(i);
    }
    return r_node_list;
  }
};

RCPP_MODULE(mod) {
  class_<Graph>("Graph")

      .constructor<int>()

      .method("connect", &Graph::connect)
      .method("cut", &Graph::cut)
      .method("getConnectedNodes", &Graph::getConnectedNodes);
}
