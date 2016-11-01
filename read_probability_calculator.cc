#include "read_probability_calculator.h"
#include <algorithm>
#include <cmath>
#include <cassert>

double SingleReadProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, ProbabilityChange& prob_change) {
  prob_change.added_paths.clear();
  prob_change.removed_paths.clear();
  prob_change.added_alignments.clear();
  prob_change.removed_alignments.clear();
  ComparePathSets(old_paths_, paths, prob_change.added_paths, prob_change.removed_paths);

  prob_change.new_paths_length = GetPathsLength(paths);
  prob_change.new_paths = paths;

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
  for (size_t i = 0; i < prob_change.removed_paths.size(); i++) {
    auto &p = prob_change.removed_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.removed_alignments.insert(prob_change.removed_alignments.end(), als.begin(), als.end()); 
    printf("\rdone %d/%d evals", (int) i+1, (int) prob_change.removed_paths.size()); 
    fflush(stdout);
  }
  printf("\n");
}

double SingleReadProbabilityCalculator::EvalTotalProbabilityFromChange(
    const ProbabilityChange& prob_change, bool write) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  // (read_id, prob_change)
  vector<pair<int, double>> changes;
  for (auto &a: prob_change.added_alignments) {
    changes.push_back(make_pair(a.read_id, GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size())));
  }
  for (auto &a: prob_change.removed_alignments) {
    changes.push_back(make_pair(a.read_id, -GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size())));
  }
  sort(changes.begin(), changes.end());
  int last_read_id = -47;
  double accumulated_prob = 0;
  for (auto &ch: changes) {
    if (ch.first != last_read_id && last_read_id != -47) {
      new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
      new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
      if (write) {
        read_probs_[last_read_id] += accumulated_prob;
      }
      accumulated_prob = 0;
    }
    accumulated_prob += ch.second;
    last_read_id = ch.first;
  }
  if (last_read_id != -47) {
    new_prob -= GetRealReadProbability(read_probs_[last_read_id], last_read_id) / read_set_->size();
    new_prob += GetRealReadProbability(read_probs_[last_read_id] + accumulated_prob, last_read_id) / read_set_->size();
    if (write) {
      read_probs_[last_read_id] += accumulated_prob;
    }
  }
  if (write) total_log_prob_ = new_prob;
  return new_prob;
}

double SingleReadProbabilityCalculator::GetAlignmentProb(
    int dist, int read_length) const {
  return pow(mismatch_prob_, dist) * pow(1 - mismatch_prob_, read_length - dist);
}

void SingleReadProbabilityCalculator::ApplyProbabilityChange(
    const ProbabilityChange& prob_change) {
  EvalTotalProbabilityFromChange(prob_change, true);
  old_paths_ = prob_change.new_paths;
  old_paths_length_ = prob_change.new_paths_length; 
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
  return max(log(max(0.0, prob)), GetMinLogProbability((*read_set_)[read_id].size()));
}

int SingleReadProbabilityCalculator::GetPathsLength(const vector<Path>& paths) const {
  int ret = 0;
  for (auto &p: paths) {
    ret += p.ToString(true).size();
  }
  return ret;
}

GlobalProbabilityCalculator::GlobalProbabilityCalculator(const Config& config) {
  for (auto &single_short_reads: config.single_short_reads()) {
    ReadSet<>* rs = new ReadSet<>();
    rs->LoadReadSet(single_short_reads.filename());
    read_sets_.push_back(rs);
    single_read_calculators_.push_back(make_pair(SingleReadProbabilityCalculator(
          rs, single_short_reads.mismatch_prob(),
          single_short_reads.min_prob_start(),
          single_short_reads.min_prob_per_base(),
          single_short_reads.penalty_constant(),
          single_short_reads.penalty_step()), single_short_reads.weight()));
  }
}

double GlobalProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, ProbabilityChanges& prob_changes) {
  prob_changes.single_read_changes.clear();
  double total_prob = 0;
  for (auto &single_read_calculator: single_read_calculators_) {
    ProbabilityChange ch;
    double prob = single_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * single_read_calculator.second;
    prob_changes.single_read_changes.push_back(ch);
  }
  return total_prob;
}

void GlobalProbabilityCalculator::ApplyProbabilityChanges(
    const ProbabilityChanges& prob_changes) {
  assert(prob_changes.single_read_changes.size() == single_read_calculators_.size());
  for (size_t i = 0; i < single_read_calculators_.size(); i++) {
    single_read_calculators_[i].first.ApplyProbabilityChange(prob_changes.single_read_changes[i]);
  }
}
