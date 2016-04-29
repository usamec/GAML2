#ifndef MOVES_H__
#define MOVES_H__

#include "path.h"

class MoveConfig {
 public:
  int big_node_threshold;
  int rand_extend_step_threshold;
  int rand_extend_distance_threshold;

  MoveConfig() :
    big_node_threshold(500),
    rand_extend_step_threshold(50),
    rand_extend_distance_threshold(1000)
    {}
};

void MakeMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config);
bool TryMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config);

#endif
