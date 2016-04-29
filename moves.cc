#include "moves.h"
#include <algorithm>
#include <cassert>

int FindPathWithSameEnding(vector<Path>& paths, int pi) {
  for (size_t i = 0; i < paths.size(); i++) {
    auto &p = paths[i];
    auto rp = p.GetReverse();
    if (i == pi) continue;
    if (p[0] == paths[pi].back()) {
      return i;
    }
    if (rp[0] == paths[pi].back()) {
      paths[i] = rp;
      return i;
    }
  }
  return -1;
}

void MergePaths(vector<Path>& paths, int pi, int same_end) {
  assert(paths[pi].back() == paths[same_end][0]);
  paths[pi].AppendPath(paths[same_end], 1);
  swap(paths[same_end], paths[paths.size()-1]);
  paths.pop_back();
}

bool ExtendPathsRandomly(const vector<Path>& paths, vector<Path>& out_paths,
                         const MoveConfig& config) {
  int pi = rand()%paths.size();
  out_paths = paths;
  random_shuffle(out_paths.begin(), out_paths.end());
  if (rand()%2 == 1) {
    out_paths[pi].Reverse();
  }
  if (!out_paths[pi].ExtendRandomly(config.big_node_threshold,
                                    config.rand_extend_step_threshold,
                                    config.rand_extend_distance_threshold)) {
    return false;
  }
  int same_end = FindPathWithSameEnding(out_paths, pi);
  if (same_end == -1) {
    return true;
  }

  MergePaths(out_paths, pi, same_end);
  return true;
}

void MakeMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config) {
  while (true) {
    out_paths.clear();
    if (TryMove(paths, out_paths, config)) return;
  }
}

bool TryMove(const vector<Path>& paths, vector<Path>& out_paths, const MoveConfig& config) {
  int move = rand()%1;
  if (move == 0) {
    return ExtendPathsRandomly(paths, out_paths, config);
  }
  return false;
}
