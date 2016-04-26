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
}

double SingleReadProbabilityCalculator::EvalTotalProbabilityFromChange(
    const ProbabilityChange& prob_change) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  // (read_id, prob_change)
  vector<pair<int, double>> changes;
  for (auto &a: prob_change.added_alignments) {
    changes.push_back(make_pair(a.read_id, GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size())));
  }
  sort(changes.begin(), changes.end());
  int last_read_id = -47;
  double accumulated_prob = 0;
  for (auto &ch: changes) {
    if (ch.first != last_read_id && last_read_id != -47) {
      new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
      new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
      accumulated_prob = 0;
    }
    accumulated_prob += ch.second;
    last_read_id = ch.first;
  }
  if (last_read_id != -47) {
    new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
    new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
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
    read_probs_[i] = 0;
    ret += GetMinLogProbability((*read_set_)[i].size()) / read_set_->size();
  }
  return ret;
}

double SingleReadProbabilityCalculator::GetMinLogProbability(int read_length) const {
  return min_prob_start_ + read_length * min_prob_per_base_; 
}

double SingleReadProbabilityCalculator::GetRealReadProbability(double prob, int read_id) const {
  return max(log(prob), GetMinLogProbability((*read_set_)[read_id].size()));
}

int SingleReadProbabilityCalculator::GetPathsLength(const vector<Path>& paths) const {
  int ret = 0;
  for (auto &p: paths) {
    ret += p.ToString(true).size();
  }
  return ret;
}
