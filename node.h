#ifndef NODE_H__
#define NODE_H__

#include <string>
#include <fstream>
#include <vector>

using namespace std;

class Node {
 public:
  Node() {}
  Node(const string& str, int id) : str_(str), id_(id), rc_(NULL) {}

  void AddNext(Node *x) {
    next_.push_back(x);
  }

  bool IsBig(size_t threshold) const {
    return str_.size() >= threshold;
  }

  string str_;
  int id_;
  Node* rc_;
  vector<Node*> next_;
};

pair<Node*,Node*> LoadNode(istream& is, int number);

#endif
