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

  Node*& back() {
    return nodes_.back();
  }

  bool CheckPath() const;

  void AppendPath(const Path& p, int p_start=0);
  void AppendPathWithGap(const Path &p, int gap_length, int p_start=0);

  void Reverse();
  Path GetReverse();

  string ToDebugString() const;

  string ToString(bool with_endings=false) const;

  bool IsSame(const Path& p) const;

  bool operator==(const Path &p) const {
    return nodes_ == p.nodes_;
  }

  bool ExtendRandomly(int big_node_threshold, int step_threshold, int distance_threshold);
};

vector<Path> BuildPathsFromSingleNodes(const vector<Node*>& nodes);

string PathsToDebugString(const vector<Path>& paths);

void PathsToFasta(const vector<Path>& paths, ostream &of);

void ComparePathSets(const vector<Path>& a,
                     const vector<Path>& b,
                     vector<Path>& added,
                     vector<Path>& removed);
#endif 
