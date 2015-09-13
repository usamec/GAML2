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

  size_t size() const {
    return nodes_.size();
  }

  bool CheckPath() const;

  void AppendPath(const Path& p);

  void Reverse();
  Path GetReverse();

  string ToDebugString() const;
};

vector<Path> BuildPathsFromSingleNodes(const vector<Node*>& nodes);

string PathsToDebugString(const vector<Path>& paths);

#endif 