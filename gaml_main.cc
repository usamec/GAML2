#include "graph.h"
#include "path.h"
#include <iostream>

int main(int argc, char** argv) {
  Graph *g = LoadGraph(argv[1]);
  cout << "Load graph with " << g->nodes_.size() << " nodes" << endl;

  int threshold = 500;

  vector<Path> paths = BuildPathsFromSingleNodes(g->GetBigNodes(threshold));

  cout << PathsToDebugString(paths) << endl;
}
