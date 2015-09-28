#ifndef PATH_ALIGNER_H__
#define PATH_ALIGNER_H__

#include "read_set.h"
#include "hash_util.h"
#include "path.h"

class PathAligner {
 public:
  PathAligner() {}
  PathAligner(ReadSet<>* read_set) : read_set_(read_set) {}

  vector<ReadAlignment> GetAlignmentsForPath(const Path& p);

  ReadSet<>* read_set_;
};

#endif
