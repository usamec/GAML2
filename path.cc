#include "path.h"
#include <algorithm>
#include <sstream>

vector<Path> BuildPathsFromSingleNodes(const vector<Node*>& nodes) {
  vector<Path> ret;
  for (auto &n: nodes) {
    ret.push_back(Path({n}));
  }
  return ret;
}

bool Path::CheckPath() const {
  for (size_t i = 0; i + 1 < nodes_.size(); i++) {
    if (find(nodes_[i]->next_.begin(), nodes_[i]->next_.end(), nodes_[i+1]) ==
        nodes_[i]->next_.end()) {
      return false;
    }    
  }
  return true;
}

string Path::ToDebugString() const {
  stringstream ret;
  ret << "(";
  for (size_t i = 0; i < nodes_.size(); i++) {
    ret << nodes_[i]->id_;
    if (i + 1 != nodes_.size()) {
      ret << ",";
    }
  }
  ret << ")";
  return ret.str();
}

string PathsToDebugString(const vector<Path>& paths) {
  stringstream ret;
  ret << paths.size() << " paths:  ";
  for (size_t i = 0; i < paths.size(); i++) {
    ret << paths[i].ToDebugString();
    if (i + 1 != paths.size()) {
      ret << "   ";
    }
  }
  return ret.str();
}
