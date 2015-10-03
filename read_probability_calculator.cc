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
  for (size_t i = 0; i < prob_change.added_paths.size(); i++) {
    auto &p = prob_change.added_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.added_alignments.insert(prob_change.added_alignments.end(), als.begin(), als.end()); 
    printf("\rdone %d/%d evals", (int) i+1, (int) prob_change.added_paths.size()); 
    fflush(stdout);
  }
  printf("\n");
  sort(prob_change.added_alignments.begin(), prob_change.added_alignments.end());
}

double SingleReadProbabilityCalculator::EvalTotalProbabilityFromChange(
    const ProbabilityChange& prob_change) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  int last_read_id = -47;
  double accumulated_prob = 0;
  for (size_t i = 0; i < prob_change.added_alignments.size(); i++) {
    auto &a = prob_change.added_alignments[i];
    if (a.read_id != last_read_id && last_read_id != -47) {
      new_prob -= read_probs_[last_read_id] / read_set_->size();
      new_prob += log(accumulated_prob) / read_set_->size();
      accumulated_prob = 0;
    }
    accumulated_prob += GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size());
    last_read_id = a.read_id;
  }
  if (last_read_id != -47) {
    new_prob -= read_probs_[last_read_id] / read_set_->size();
    new_prob += log(accumulated_prob) / read_set_->size();
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
    ret += read_probs_[i] / read_set_->size();
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
