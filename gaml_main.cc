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
#include <csignal>

bool NOT_SIGINTED = true;
void signalHandler(int signum) {
  NOT_SIGINTED = false;
}

void PerformOptimization(GlobalProbabilityCalculator& probability_calculator,
                         const Config& gaml_config, vector<Path>& paths) {
  // @TODO add random seed to config? (i14)
  const int random_seed = 47;
  cout << "random seed: " << random_seed << endl;
  default_random_engine generator(random_seed);

  ProbabilityChanges prob_changes;
  double old_prob = probability_calculator.GetPathsProbability(paths, prob_changes);
  cout << "starting probability: " << old_prob << endl;

  probability_calculator.CommitProbabilityChanges(prob_changes);

  cout << "INITIAL PATHS:" << endl;
  cout << PathsToDebugString(paths) << endl;

  MoveConfig move_config;
  for (int it_num = 1; it_num <= gaml_config.num_iterations() && NOT_SIGINTED; it_num++) {
    double T = gaml_config.t0() / log(it_num / gaml_config.n_divisor() + 1);
    int total_size = 0;
    for (auto &p: paths) total_size += p.ToString(true).size();
    cout << "ITERATION: " << it_num << "\tT: " << T << "\tPROB: " << old_prob << "\t SIZE: " << total_size <<  endl;

    vector<Path> new_paths;
    bool accept_high_prob;
    MakeMove(paths, new_paths, move_config, probability_calculator, accept_high_prob);
    double new_prob = probability_calculator.GetPathsProbability(new_paths, prob_changes);
    cout << "PROPOSED PROB: " << new_prob;

    bool accept = false;
    if (new_prob > old_prob) {
      accept = true;
      cout << " better than old; ";
    } else if (accept_high_prob) {
      double prob = exp((new_prob - old_prob) / T);
      uniform_real_distribution<double> dist(0.0, 1.0);
      double samp = dist(generator);
      if (samp < prob) {
        cout << " worse, but temperature; ";
        accept = true;
      }
    }

    if (accept) {
      cout << "\tACCEPT";
      if (old_prob < new_prob){
        // debug
        for (auto &np: new_paths) {
          bool is_new = true;
          for (auto &p: paths) {
            if (p.IsSame(np)) {
              is_new = false;
              break;
            }
          }
          if (is_new) np.history_ += PATH_BETTER;
        }
      }
      old_prob = new_prob;

      paths = new_paths;
      probability_calculator.CommitProbabilityChanges(prob_changes);
    }
    cout << endl << "PROPOSED PATHS:\n" <<  PathsToDebugString(new_paths) << endl << endl;

    // continual output
    if (it_num % gaml_config.output_flush_freq() == 0) {
      ofstream of(gaml_config.output_file());
      PathsToFasta(paths, of);
    }
  }

  // final output
  ofstream of(gaml_config.output_file());
  PathsToFasta(paths, of);
}

void global_run_logging(const string& log_filename, int argc, char** argv) {
  fstream fs;
  fs.open(log_filename, fstream::out | fstream::app);
  fs << "Running GAML2..." << endl;
  const auto ctt = time(0);
  fs << "TIMESTAMP: " << asctime(localtime(&ctt)); // contains \n symbol
  fs << "COMMAND: ";
  for (int i = 0; i < argc; i++) {
    fs << argv[i] << " ";
  }
  fs << endl;

  fstream cfg_file;
  cfg_file.open(argv[1], fstream::in);
  string buff = "x";
  while (cfg_file >> buff) {
    fs << "\t" << buff << endl;
  }

  cfg_file.close();
  fs.close();
}

int main(int argc, char** argv) {
  // Adding correct handling of Ctrl+C signal (SIGINT)
  signal(SIGINT, signalHandler);

  ifstream config_file(argv[1]);
  google::protobuf::io::IstreamInputStream config_stream(&config_file);

  Config gaml_config;
  google::protobuf::TextFormat::Parse(&config_stream, &gaml_config);

  global_run_logging(gaml_config.global_run_log_file(), argc, argv);

  cout << gaml_config.starting_graph() << endl; 

  Graph *g = LoadGraph(gaml_config.starting_graph());
  cout << "Loaded graph with " << g->nodes_.size() << " nodes" << endl;

  GlobalProbabilityCalculator probability_calculator(gaml_config);

  vector<Path> paths = BuildPathsFromSingleNodes(g->GetBigNodes(gaml_config.big_nodes_threshold()));

  cout << PathsToDebugString(paths) << endl;

  PerformOptimization(probability_calculator, gaml_config, paths);
}
