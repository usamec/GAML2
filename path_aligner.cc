#include "path_aligner.h"

vector<SingleReadAlignment> SingleShortReadPathAligner::GetAlignmentsForPath(const Path& p) {
  vector<SingleReadAlignment> ret;

  string genome = p.ToString(true);

  ret = single_short_read_set_->GetAlignments(genome);
  return ret;
}
