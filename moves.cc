#include "moves.h"

void MakeMove(const vector<Path>& paths, vector<Path>& out_paths) {
  while (true) {
    out_paths.clear();
    if (TryMove(paths, out_paths)) return;
  }
}

bool TryMove(const vector<Path>& paths, vector<Path>& out_paths) {
  return true;
}

