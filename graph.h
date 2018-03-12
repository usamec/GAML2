#ifndef GRAPH_H__
#define GRAPH_H__

#include "node.h"
#include <tuple>

class Graph {
 public:
  Graph() {}

  vector<Node*> GetBigNodes(int threshold) const;

  vector<Node*> ReachForwardWithThreshold(Node* start, int threshold) const;
  vector<Node*> ReachLocalWithThreshold(Node* start, int threshold) const;
  vector<Node*> GetDrainageBasinForNode(Node *target);

  vector<Node*> nodes_;
  // `k_` := k from de Bruijn graph (velvet output graph header, see `::Load Graph` )
  int k_;
};

int VelvetNumToMyNum(int x);
pair<int, int> LoadArc(istream &is);
Graph* LoadGraph(const string& filename);
Graph* LoadGraph(istream &is);

#endif
