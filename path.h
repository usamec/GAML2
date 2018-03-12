#ifndef PATH_H__
#define PATH_H__

#include "node.h"

static const string PATH_BIGNODE = "I";
static const string PATH_EXTEND_RANDOMLY = "E";
static const string PATH_BREAK = "B";
static const string PATH_JOIN = "J";
static const string PATH_BETTER = "+";
static const string PATH_UNTANGLE = "U";

class Path {
 public:
  Path() {}
  explicit Path(const vector<Node*> nodes) : nodes_(nodes) {}
  explicit Path(const vector<Node*> nodes, const string history) : nodes_(nodes), history_(history) {}
  vector<Node*> nodes_;
  // for debug purposes
  string history_;

  Node*& operator[](size_t x) {
    return nodes_[x];
  }

  Path& operator=(const Path& other) {
    nodes_ = other.nodes_;
    history_ = other.history_;
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
  Path GetReverse() const;

  string ToDebugString() const;

  string ToString(bool with_endings=false) const;

  bool IsSame(const Path& p) const;

  bool IsSameNoReverse(const Path& p) const;

  bool isDisjoint(const Path& p) const;

  // check for disjoin but only for ends
  bool isPartlyDisjoint(const Path &p, int dist) const;

  bool operator==(const Path &p) const {
    return nodes_ == p.nodes_;
  }

  bool ExtendRandomly(int big_node_threshold, int step_threshold, int distance_threshold);

  // Split path into <0, pos) and <pos, ...) and removes small nodes from ends
  Path CutAt(int pos, int big_node_threshold);
};

vector<Path> BuildPathsFromSingleNodes(const vector<Node*>& nodes);

string PathsToDebugString(const vector<Path>& paths);

void PathsToFasta(const vector<Path>& paths, ostream &of);

void ComparePathSets(const vector<Path>& a,
                     const vector<Path>& b,
                     vector<Path>& added,
                     vector<Path>& removed);
#endif 
