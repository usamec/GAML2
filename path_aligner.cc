#include "path_aligner.h"

vector<ReadAlignment> PathAligner::GetAlignmentsForPath(const Path& p) {
  vector<ReadAlignment> ret;

  string genome = p.ToString(true);

  ret = read_set_->GetAlignments(genome);
  return ret;
}
