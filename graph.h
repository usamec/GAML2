#ifndef GRAPH_H__
#define GRAPH_H__

#include "node.h"
#include <tuple>

class Graph {
 public:
  Graph() {}

  vector<Node*> nodes_;
  int k_;
};

int VelvetNumToMyNum(int x);
pair<int, int> LoadArc(istream &is);
Graph* LoadGraph(string& filename);
Graph* LoadGraph(istream &is);

#endif
