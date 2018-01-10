#include "path_aligner.h"

vector<SingleReadAlignment> SingleShortReadPathAligner::GetAlignmentsForPath(const Path& p) {
  vector<SingleReadAlignment> ret;

  string genome = p.ToString(true);

  ret = single_short_read_set_->GetAlignments(genome);
  return ret;
}
vector<PairedReadAlignment> PairedReadPathAligner::GetAlignmentsForPath(const Path &p) {
  vector<PairedReadAlignment> ret;
  string genome = p.ToString(true);
  ret = paired_read_set_->GetAlignments(genome);
  return ret;
}
