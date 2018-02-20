#include "path_aligner.h"

vector<SingleReadAlignment> SingleShortReadPathAligner::GetAlignmentsForPath(const Path& p) {
  vector<SingleReadAlignment> ret;

  const int pos = GetAlignmentPos(p);
  if (pos == -1) {
    string genome = p.ToString(true);
    ret = single_short_read_set_->GetAlignments(genome);
    InsertAlignmentForPath(p, ret);
  }
  else {
    ret = GetCachedAlignmentByPos(pos);
  }

  return ret;
}
vector<PairedReadAlignment> PairedReadPathAligner::GetAlignmentsForPath(const Path &p) {
  vector<PairedReadAlignment> ret;

  const int pos = GetAlignmentPos(p);
  if (pos == -1) {
    const vector<SingleReadAlignment> als1 = GetPartAlignmentsForPath(p, 0);
    const vector<SingleReadAlignment> als2 = GetPartAlignmentsForPath(p, 1);

    // assuming als1 and als2 are sorted by read position as primary key

    auto it1 = als1.begin();
    auto it2 = als2.begin();

    const auto &reads_1_ = left_aligner_.single_short_read_set_;
    const auto &reads_2_ = right_aligner_.single_short_read_set_;
    const int reads_1_size = reads_1_->size();

    // run through both vectors, find groups with equal read_id and do cross-check
    // for potential pairing (with respect to the assumed orientation)
    vector<SingleReadAlignment> current_als1, current_als2;
    for (int current_read_id = 0; current_read_id < reads_1_size && it1 != als1.end() && it2 != als2.end(); current_read_id++) {
      current_als1.clear();
      current_als2.clear();
      while (it1 != als1.end() && it1->read_id < current_read_id) it1++;
      while (it1 != als1.end() && it1->read_id == current_read_id) {
        current_als1.push_back(*it1);
        it1++;
      }
      while (it2 != als2.end() && it2->read_id < current_read_id) it2++;
      while (it2 != als2.end() && it2->read_id == current_read_id) {
        current_als2.push_back(*it2);
        it2++;
      }
      //printf("c12 %d,%d\t", (int)current_als1.size(), (int)current_als2.size());
      if (current_als1.size() > 0 && current_als2.size() > 0) {
        // cross check for finding paired alignments
        string orient;
        int insert_length;
        for (auto &a1: current_als1) {
          for (auto &a2: current_als2) {
            tie(orient, insert_length) = eval_orientation(a1, reads_1_[current_read_id].size(), a2, reads_2_[current_read_id].size());
            ret.push_back(PairedReadAlignment(a1, a2, orient, insert_length));
          }
        }
      }
    }

    InsertAlignmentForPath(p, ret);
  }
  else {
    ret = GetCachedAlignmentByPos(pos);
  }

  return ret;
}
vector<SingleReadAlignment> PairedReadPathAligner::GetPartAlignmentsForPath(const Path &p, int part) {
  vector<SingleReadAlignment> ret;
  if (part == 0) {
    ret = left_aligner_.GetAlignmentsForPath(p);
  }
  else if (part == 1) {
    ret = right_aligner_.GetAlignmentsForPath(p);
  }
  return ret;
}
