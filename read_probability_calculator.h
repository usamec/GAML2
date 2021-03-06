#ifndef READ_PROBABILITY_CALCULATOR_H__
#define READ_PROBABILITY_CALCULATOR_H__

#include "path_aligner.h"
#include "config.pb.h"

struct ProbabilityChange {
  vector<Path> added_paths;
  vector<Path> removed_paths;

  vector<ReadAlignment> added_alignments;
  vector<ReadAlignment> removed_alignments;

  int new_paths_length;

  vector<Path> new_paths;
};

struct ProbabilityChanges {
  vector<ProbabilityChange> single_read_changes;
};

class SingleReadProbabilityCalculator {
 public:
  SingleReadProbabilityCalculator(
      ReadSet<>* read_set, double mismatch_prob,
      double min_prob_start, double min_prob_per_base,
      double penalty_constant, int penalty_step) :
        read_set_(read_set), path_aligner_(read_set),
        mismatch_prob_(mismatch_prob),
        min_prob_start_(min_prob_start), min_prob_per_base_(min_prob_per_base),
        penalty_constant_(penalty_constant), penalty_step_(penalty_step),
        old_paths_length_(1) {
    read_probs_.resize(read_set_->size());
    total_log_prob_ = InitTotalLogProb();
  }

  // Call this first
  double GetPathsProbability(
      const vector<Path>& paths, ProbabilityChange& prob_change);

  // Call this after you are happy with current result (i.e. you got better
  // probability)
  void ApplyProbabilityChange(const ProbabilityChange& prob_change);

 private:
  double InitTotalLogProb();

  double GetMinLogProbability(int read_length) const;

  // max(min_prob, prob)
  double GetRealReadProbability(double prob, int read_id) const;

  // Evals change with filled added and removed paths
  void EvalProbabilityChange(ProbabilityChange& prob_change);

  // Get total probability from change and cached data
  double EvalTotalProbabilityFromChange(const ProbabilityChange& prob_change, bool write=false);

  int GetPathsLength(const vector<Path>& paths) const;

  double GetAlignmentProb(int dist, int read_length) const;

  ReadSet<>* read_set_;
  PathAligner path_aligner_;
  double mismatch_prob_;
  double min_prob_start_;
  double min_prob_per_base_;
  double penalty_constant_;
  int penalty_step_;

  vector<double> read_probs_;
  double total_log_prob_;
  int old_paths_length_;
  vector<Path> old_paths_;
};

class GlobalProbabilityCalculator {
 public:
  GlobalProbabilityCalculator(const Config &config);

  // Call this first
  double GetPathsProbability(
      const vector<Path>& paths, ProbabilityChanges& prob_changes);

  // Call this after you are happy with current result (i.e. you got better
  // probability)
  void ApplyProbabilityChanges(const ProbabilityChanges& prob_changes);

 private:
  vector<ReadSet<>*> read_sets_;
  // (prob calculator, weight)
  vector<pair<SingleReadProbabilityCalculator, double>> single_read_calculators_;

};

#endif
