#ifndef PATH_H__
#define PATH_H__

#include "node.h"

class Path {
 public:
  Path() {}
  Path(const vector<Node*> nodes) : nodes_(nodes) {}
  vector<Node*> nodes_;

  Node*& operator[](size_t x) {
    return nodes_[x];
  }

  bool CheckPath() const;
};

vector<Path> BuildPathsFromSingleNodes(const vector<Node*>& nodes);

#endif 
