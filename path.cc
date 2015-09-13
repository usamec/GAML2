#include "path.h"
#include <algorithm>

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
