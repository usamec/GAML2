#include "graph.h"

int main(int argc, char** argv) {
  Graph *g = LoadGraph(argv[1]);
  int threshold = 500;
  int start = atoi(argv[2]);
  vector<Node*> nodes = g->ReachLocalWithThreshold(g->nodes_[start], threshold);
  printf("digraph G {\n");
  for (auto &n: nodes) {
    for (auto &next: n->next_) {
      printf("\"%d(%lu)\" -> \"%d(%lu)\"\n", n->id_, n->str_.size(), next->id_, next->str_.size());
    }
    for (auto &nextr: n->rc_->next_) {
      auto &next = nextr->rc_;
      printf("\"%d(%lu)\" -> \"%d(%lu)\"\n", next->id_, next->str_.size(), n->id_, n->str_.size());
    }
  }
  printf("}\n");
} 
