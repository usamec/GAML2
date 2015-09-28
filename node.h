#ifndef NODE_H__
#define NODE_H__

#include <string>
#include <fstream>
#include <vector>

using namespace std;

class Graph;

class Node {
 public:
  Node() {}
  Node(const string& str, int id, Graph* graph) :
      str_(str), id_(id), rc_(NULL), graph_(graph) {}

  void AddNext(Node *x) {
    next_.push_back(x);
  }

  bool IsBig(size_t threshold) const {
    return str_.size() >= threshold;
  }

  bool IsGap() const {
    return id_ < 0;
  }

  int GapLength() const {
    return id_ < 0 ? -id_ : 0;
  }

  string str_;
  int id_;
  Node* rc_;
  vector<Node*> next_;
  Graph* graph_;
};

pair<Node*,Node*> LoadNode(istream& is, int number, Graph* gr);
Node* MakeGap(int gap_length);

#endif
