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
#include <cmath>
#include <random>

void PerformOptimization(GlobalProbabilityCalculator& probability_calculator,
                         const Config& gaml_config, vector<Path>& paths) {
  // @TODO add random seed to config? (i14)
  default_random_engine generator(47);
  ProbabilityChanges prob_changes;
  // SO FAR SO GOOD
  double old_prob = probability_calculator.GetPathsProbability(paths, prob_changes);
  cout << "starting probability: " << old_prob << endl;
  probability_calculator.CommitProbabilityChanges(prob_changes);

  cout << PathsToDebugString(paths) << endl;
  MoveConfig move_config;
  for (int it_num = 1; it_num <= gaml_config.num_iterations(); it_num++) {
    double T = gaml_config.t0() / log(it_num / gaml_config.n_divisor() + 1);
    cout << "Iter: " << it_num << " T: " << T << endl;

    vector<Path> new_paths;
    bool accept_high_prob;
    MakeMove(paths, new_paths, move_config, accept_high_prob);
    double new_prob = probability_calculator.GetPathsProbability(new_paths, prob_changes);
    cout << new_prob;

    bool accept = false;
    if (new_prob > old_prob) {
      accept = true;
    } else if (accept_high_prob) {
      double prob = exp((new_prob - old_prob) / T);
      uniform_real_distribution<double> dist(0.0, 1.0);
      double samp = dist(generator);
      if (samp < prob) {
        cout << " higher";
        accept = true;
      }
    }

    if (accept) {
      cout << " accept";
      old_prob = new_prob;
      paths = new_paths;
      probability_calculator.CommitProbabilityChanges(prob_changes);
    }
    cout << endl << PathsToDebugString(new_paths) << endl << endl;
  }

  ofstream of(gaml_config.output_file());
  PathsToFasta(paths, of);
}

int main(int argc, char** argv) {
  ifstream config_file(argv[1]);
  google::protobuf::io::IstreamInputStream config_stream(&config_file);

  Config gaml_config;
  google::protobuf::TextFormat::Parse(&config_stream, &gaml_config);

  cout << gaml_config.starting_graph() << endl; 

  Graph *g = LoadGraph(gaml_config.starting_graph());
  cout << "Loaded graph with " << g->nodes_.size() << " nodes" << endl;

  GlobalProbabilityCalculator probability_calculator(gaml_config);

  vector<Path> paths = BuildPathsFromSingleNodes(g->GetBigNodes(gaml_config.big_nodes_threshold()));

  cout << PathsToDebugString(paths) << endl;

  PerformOptimization(probability_calculator, gaml_config, paths);
}
