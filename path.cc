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
    if (nodes_[i]->IsGap() || nodes_[i+1]->IsGap()) continue;
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

void Path::AppendPath(const Path& p) {
  nodes_.insert(nodes_.end(), p.nodes_.begin(), p.nodes_.end());
}

void Path::AppendPathWithGap(const Path& p, int gap_length) {
  nodes_.push_back(MakeGap(gap_length));
  nodes_.insert(nodes_.end(), p.nodes_.begin(), p.nodes_.end());
}

void Path::Reverse() {
  vector<Node*> nodes_new;
  for (int i = nodes_.size() - 1; i >= 0; i--) {
    nodes_new.push_back(nodes_[i]->rc_);
  }
  nodes_ = nodes_new;
}

Path Path::GetReverse() {
  vector<Node*> nodes_new;
  for (int i = nodes_.size() - 1; i >= 0; i--) {
    nodes_new.push_back(nodes_[i]->rc_);
  }
  Path ret;
  ret.nodes_ = nodes_new;
  return ret;
}
