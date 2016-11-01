#include "graph.h"
#include "path.h"
#include "read_set.h"
#include "read_probability_calculator.h"
#include "moves.h"
#include "config.pb.h"
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
  ifstream config_file(argv[1]);
  google::protobuf::io::IstreamInputStream config_stream(&config_file);

  Config gaml_config;
  google::protobuf::TextFormat::Parse(&config_stream, &gaml_config);

  cout << gaml_config.starting_graph() << endl; 

  Graph *g = LoadGraph(gaml_config.starting_graph());
  cout << "Loaded graph with " << g->nodes_.size() << " nodes" << endl;

  GlobalProbabilityCalculator probability_calculator(gaml_config);

  int threshold = 500;

  vector<Path> paths = BuildPathsFromSingleNodes(g->GetBigNodes(threshold));

  cout << PathsToDebugString(paths) << endl;

  ProbabilityChanges prob_changes;
  double old_prob = probability_calculator.GetPathsProbability(paths, prob_changes);
  cout << old_prob << endl;
  probability_calculator.ApplyProbabilityChanges(prob_changes);

  ofstream of(argv[3]);
  PathsToFasta(paths, of);

  cout << PathsToDebugString(paths) << endl;
  MoveConfig config;
  for (int i = 0; i < 10; i++) {
    cout << "Iter: " << i << endl;
    vector<Path> new_paths;
    MakeMove(paths, new_paths, config);
    double new_prob = probability_calculator.GetPathsProbability(new_paths, prob_changes);
    cout << new_prob;
    if (new_prob > old_prob) {
      cout << " accept";
      old_prob = new_prob;
      paths = new_paths;
      probability_calculator.ApplyProbabilityChanges(prob_changes);
    }
    cout << endl << PathsToDebugString(new_paths) << endl << endl;
  }
}
