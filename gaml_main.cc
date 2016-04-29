#include "graph.h"
#include "path.h"
#include "read_set.h"
#include "read_probability_calculator.h"
#include "moves.h"
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
  double old_prob = pc.GetPathsProbability(paths, prob_change);
  cout << old_prob << endl;
//  pc.ApplyProbabilityChange(prob_change);

  ofstream of(argv[3]);
  PathsToFasta(paths, of);

  cout << PathsToDebugString(paths) << endl;
  MoveConfig config;
  for (int i = 0; i < 10; i++) {
    cout << "Iter: " << i << endl;
    vector<Path> new_paths;
    MakeMove(paths, new_paths, config);
    double new_prob = pc.GetPathsProbability(new_paths, prob_change);
    cout << new_prob;
    if (new_prob > old_prob) {
      cout << " accept";
      old_prob = new_prob;
      paths = new_paths;
//      pc.ApplyProbabilityChange(prob_change);
    }
    cout << endl << PathsToDebugString(new_paths) << endl << endl;
  }
}
