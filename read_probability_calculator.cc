#include "read_probability_calculator.h"
#include <algorithm>
#include <cmath>
#include <cassert>

double SingleReadProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, SingleProbabilityChange& prob_change) {
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
    SingleProbabilityChange& prob_change) {
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
    const SingleProbabilityChange& prob_change, bool write) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  // (read_id, prob_change)
  vector<pair<int, double>> changes;
  for (auto &a: prob_change.added_alignments) {
    //changes.push_back(make_pair(a.read_id, GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size())));
    changes.emplace_back(a.read_id, GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size()));
  }
  for (auto &a: prob_change.removed_alignments) {
    //changes.push_back(make_pair(a.read_id, -GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size())));
    changes.emplace_back(a.read_id, -GetAlignmentProb(a.dist, (*read_set_)[a.read_id].size()));
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

inline double SimpleReadProbModel (int dist, int read_length, double mismatch_prob) {
  // @DONE corrected epsilon for e/3 (error from the GAML article)
  return pow(mismatch_prob/3, dist) * pow(1 - mismatch_prob, read_length - dist);
}

double SingleReadProbabilityCalculator::GetAlignmentProb(
    int dist, int read_length) const {
  return SimpleReadProbModel(dist, read_length, mismatch_prob_);
}

void SingleReadProbabilityCalculator::CommitProbabilityChange(
    const SingleProbabilityChange &prob_change) {
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

double PairedReadProbabilityCalculator::InitTotalLogProb() {
  double ret = 0;
  for (size_t i = 0; i < read_set_->size(); i++) {
    read_probs_[i] = 0;
    ret += GetMinLogProbability((*read_set_)[i].first.size() + (*read_set_)[i].second.size()) / read_set_->size();
  }
  return ret;
}
void PairedReadProbabilityCalculator::CommitProbabilityChange(const PairedProbabilityChange &prob_change) {
  EvalTotalProbabilityFromChange(prob_change, true);
  old_paths_ = prob_change.new_paths;
  old_paths_length_ = prob_change.new_paths_length;

}
double PairedReadProbabilityCalculator::GetPathsProbability(const vector<Path> &paths,
                                                            PairedProbabilityChange &prob_change) {
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
void PairedReadProbabilityCalculator::EvalProbabilityChange(PairedProbabilityChange &prob_change, bool debug_output) {
  for (size_t i = 0; i < prob_change.added_paths.size(); i++) {
    auto &p = prob_change.added_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.added_alignments.insert(prob_change.added_alignments.end(), als.begin(), als.end());
    if (debug_output) printf("\n(added_paths) done %d/%d evals; %d alignments", (int) i+1, (int) prob_change.added_paths.size(), (int) als.size());
    if (debug_output) fflush(stdout);
  }
  if (debug_output) printf("\n");
  for (size_t i = 0; i < prob_change.removed_paths.size(); i++) {
    auto &p = prob_change.removed_paths[i];
    auto als = path_aligner_.GetAlignmentsForPath(p);
    prob_change.removed_alignments.insert(prob_change.removed_alignments.end(), als.begin(), als.end());
    if (debug_output) printf("\n(removed_paths) done %d/%d evals;  %d alignments", (int) i+1, (int) prob_change.removed_paths.size(), (int) als.size());
    if (debug_output) fflush(stdout);
  }
  if (debug_output) printf("\n");
}

double PairedReadProbabilityCalculator::GetMinLogProbability(int read_length) const {
  // see "Reads that have no good alignment to A" in the  GAML article
  return min_prob_start_ + read_length * min_prob_per_base_;
}

double PairedReadProbabilityCalculator::GetRealReadProbability(double prob, int read_id) const {
  return max(log(max(0.0, prob)), GetMinLogProbability((*read_set_)[read_id].first.size() + (*read_set_)[read_id].second.size()));
}

double PairedReadProbabilityCalculator::EvalTotalProbabilityFromChange(const PairedProbabilityChange &prob_change,
                                                                       bool write) {
  double new_prob = total_log_prob_;
  new_prob += log(old_paths_length_);
  new_prob -= log(prob_change.new_paths_length);

  // (read_id, prob_change)
  vector< pair<int, double> > changes;
  for (auto &a: prob_change.added_alignments) {
    //changes.push_back(make_pair(a.read_id, GetAlignmentProb(a)));
    changes.emplace_back(a.read_id, GetAlignmentProb(a));
  }

  for (auto &a: prob_change.removed_alignments) {
    //changes.push_back(make_pair(a.read_id, -GetAlignmentProb(a)));
    changes.emplace_back(a.read_id, -GetAlignmentProb(a));
  }

  sort(changes.begin(), changes.end());

  const int BLANK = -47;
  int last_read_id = BLANK;
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

int PairedReadProbabilityCalculator::GetPathsLength(const vector<Path> &paths) const {
  // same as for SingleReadProbabilityCalculator
  // @TODO refactor later
  int ret = 0;
  for (auto &p: paths) {
    ret += p.ToString(true).size();
  }
  return ret;
}

double pdf_normal(const double x, const double mean, const double std) {
  static const double inv_sqrt_2pi = 0.3989422804014327;
  const double a = (mean - x)/std;
  return inv_sqrt_2pi / std * exp(-0.5 * a * a);
}

bool is_same_orientation (const string& o1, const string& o2) {
  if (o1 == o2) {
    return true;
  }
  else if ((o1 == "FF" || o1 == "RR") && (o2 == "FF" || o2 == "RR") ) {
    return true;
  }
  else {
    return false;
  }
}

double PairedReadProbabilityCalculator::GetAlignmentProb(const PairedReadAlignment& al) const {
  if (! is_same_orientation(al.orientation, read_set_->orientation_)) {return 0;}
  const double pp1 = SimpleReadProbModel(al.al1.dist, (*read_set_)[al.read_id].first.size(), mismatch_prob_);
  const double pp2 = SimpleReadProbModel(al.al2.dist, (*read_set_)[al.read_id].second.size(), mismatch_prob_);

  return pp1 * pp2 * pdf_normal(al.insert_length, mean_distance_, std_distance_);
}

GlobalProbabilityCalculator::GlobalProbabilityCalculator(const Config& config) {
  for (auto &single_short_reads: config.single_short_reads()) {
    SingleShortReadSet<>* rs = new SingleShortReadSet<>();
    rs->LoadReadSet(single_short_reads.filename());
    single_short_read_sets_.push_back(rs);
    single_read_calculators_.push_back(make_pair(SingleReadProbabilityCalculator(
          rs, single_short_reads.mismatch_prob(),
          single_short_reads.min_prob_start(),
          single_short_reads.min_prob_per_base(),
          single_short_reads.penalty_constant(),
          single_short_reads.penalty_step()), single_short_reads.weight()));
  }


  for (auto &paired_reads: config.paired_reads()) {
    ShortPairedReadSet<>* rs = new ShortPairedReadSet<>();
    rs->LoadReadSet(paired_reads.filename1(), paired_reads.filename2(), paired_reads.orientation());
    paired_read_sets_.push_back(rs);
    paired_read_calculators_.push_back(make_pair(
        PairedReadProbabilityCalculator(
            rs,
            paired_reads.mismatch_prob(),
            paired_reads.min_prob_start(),
            paired_reads.min_prob_per_base(),
            paired_reads.penalty_constant(),
            paired_reads.penalty_step(),
            paired_reads.mean_distance(),
            paired_reads.std_distance(),
            paired_reads.use_as_advice()
        ),
        paired_reads.weight()
    ));
  }
}

double GlobalProbabilityCalculator::GetPathsProbability(
    const vector<Path>& paths, ProbabilityChanges& prob_changes) {
  double total_prob = 0;
  prob_changes.single_read_changes.clear();
  for (auto &single_read_calculator: single_read_calculators_) {
    SingleProbabilityChange ch;
    double prob = single_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * single_read_calculator.second;
    prob_changes.single_read_changes.push_back(ch);
  }

  prob_changes.paired_read_changes.clear();
  for (auto &paired_read_calculator: paired_read_calculators_) {
    PairedProbabilityChange ch;
    double prob = paired_read_calculator.first.GetPathsProbability(paths, ch);
    total_prob += prob * paired_read_calculator.second;
    prob_changes.paired_read_changes.push_back(ch);
  }

  return total_prob;
}

void GlobalProbabilityCalculator::CommitProbabilityChanges(
    const ProbabilityChanges &prob_changes) {
  assert(prob_changes.single_read_changes.size() == single_read_calculators_.size());
  for (size_t i = 0; i < single_read_calculators_.size(); i++) {
    single_read_calculators_[i].first.CommitProbabilityChange(prob_changes.single_read_changes[i]);
  }

  assert(prob_changes.paired_read_changes.size() == paired_read_calculators_.size());
  for (size_t i = 0; i < paired_read_calculators_.size(); i++) {
    paired_read_calculators_[i].first.CommitProbabilityChange(prob_changes.paired_read_changes[i]);
  }
}
