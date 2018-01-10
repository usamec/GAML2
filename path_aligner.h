#ifndef PATH_ALIGNER_H__
#define PATH_ALIGNER_H__

#include "read_set.h"
#include "hash_util.h"
#include "path.h"

// Handles caching of alignments
// @TODO implement caching (i8)

class SingleShortReadPathAligner {
 public:
  SingleShortReadPathAligner() {}
  SingleShortReadPathAligner(SingleShortReadSet<>* single_short_read_set) : single_short_read_set_(single_short_read_set) {}


  vector<SingleReadAlignment> GetAlignmentsForPath(const Path& p);

  SingleShortReadSet<>* single_short_read_set_;
};

class PairedReadPathAligner {
 public:
  PairedReadPathAligner() {}
  PairedReadPathAligner(ShortPairedReadSet<>* paired_read_set): paired_read_set_(paired_read_set) {}

  vector<PairedReadAlignment> GetAlignmentsForPath(const Path& p);

  ShortPairedReadSet<>* paired_read_set_;
};

#endif
