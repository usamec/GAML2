#include "read_probability_calculator.h"
#include <algorithm>
#include <cmath>

double SingleReadProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, ProbabilityChange& prob_change) {
  prob_change.added_paths.clear();
  prob_change.added_paths = paths;

  prob_change.new_paths_length = GetPathsLength(paths);

  EvalProbabilityChange(prob_change);

  return EvalTotalProbabilityFromChange(prob_change);
}

void SingleReadProbabilityCalculator::EvalProbabilityChange(
    ProbabilityChange& prob_change) {
  for (auto &p: prob_change.added_paths) {
    prob_change.added_alignments = path_aligner_.GetAlignmentsForPath(p); 
  }
  sort(prob_change.added_alignments.begin(), prob_change.added_alignments.end());
}

double SingleReadProbabilityCalculator::EvalTotalProbabilityFromChange(
    const ProbabilityChange& prob_change) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  int last_read_id = -47;
  double accumulated_prob = 0;
  for (auto &a: prob_change.added_alignments) {
    if (a.read_id != last_read_id && last_read_id != -47) {
      new_prob -= read_probs_[last_read_id];
      new_prob += log(accumulated_prob);
      accumulated_prob = 0;
    }
    accumulated_prob += GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size());
    last_read_id = a.read_id;
  }
  if (last_read_id != -47) {
    new_prob -= read_probs_[last_read_id];
    new_prob += log(accumulated_prob);
  }
  return new_prob;
}

double SingleReadProbabilityCalculator::GetAlignmentProb(
    int dist, int read_length) const {
  return pow(mismatch_prob_, dist) * pow(1 - mismatch_prob_, read_length - dist);
}

void SingleReadProbabilityCalculator::ApplyProbabilityChange(
    const ProbabilityChange& prob_change) {
}

double SingleReadProbabilityCalculator::InitTotalLogProb() {
  double ret = 0;
  for (size_t i = 0; i < read_set_->size(); i++) {
    read_probs_[i] = GetMinLogProbability((*read_set_)[i].size());
    ret += read_probs_[i];
  }
  return ret;
}

double SingleReadProbabilityCalculator::GetMinLogProbability(int read_length) const {
  return min_prob_start_ + read_length * min_prob_per_base_; 
}

int SingleReadProbabilityCalculator::GetPathsLength(const vector<Path>& paths) const {
  int ret = 0;
  for (auto &p: paths) {
    ret += p.ToString(true).size();
  }
  return ret;
}
