#include "graph.h"
#include <iostream>

int main(int argc, char** argv) {
  Graph *g = LoadGraph(argv[1]);
  cout << "Load graph with " << g->nodes_.size() << " nodes" << endl;
}
