#include "node.h"

pair<Node*, Node*> LoadNode(istream& is, int number, Graph* gr) {
  string header, forward, backward;
  getline(is, header);
  getline(is, forward);
  getline(is, backward);

  Node *f = new Node(forward, number*2, gr);
  Node *b = new Node(backward, number*2+1, gr);
  f->rc_ = b;
  b->rc_ = f;

  return make_pair(f, b);
}

Node* MakeGap(int gap_length) {
  Node* ret = new Node;
  ret->id_ = -gap_length;
  ret->rc_ = ret;
  return ret;
}
