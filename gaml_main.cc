#include "graph.h"
#include "path.h"
#include "read_set.h"
#include "read_probability_calculator.h"
#include <iostream>

int main(int argc, char** argv) {
  Graph *g = LoadGraph(argv[1]);
  cout << "Load graph with " << g->nodes_.size() << " nodes" << endl;

  int threshold = 500;

  vector<Path> paths = BuildPathsFromSingleNodes(g->GetBigNodes(threshold));

  cout << PathsToDebugString(paths) << endl;

  ReadSet<> rs;
  rs.LoadReadSet(argv[2]);

  cout << "read set loaded " << rs.size() << endl;

  SingleReadProbabilityCalculator pc(&rs, 0.01, -10, -0.7, 0, 0);
  ProbabilityChange prob_change;
  cout << pc.GetPathsProbability(paths, prob_change) << endl;

  ofstream of(argv[3]);
  PathsToFasta(paths, of);
}
