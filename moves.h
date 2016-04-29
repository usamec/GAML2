#ifndef MOVES_H__
#define MOVES_H__

#include "path.h"

void MakeMove(const vector<Path>& paths, vector<Path>& out_paths);
bool TryMove(const vector<Path>& paths, vector<Path>& out_paths);

#endif
