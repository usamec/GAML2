#include "path.h"
#include "graph.h"
#include "util.h"
#include <cassert>
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

void Path::AppendPath(const Path& p, int p_start) {
  nodes_.insert(nodes_.end(), p.nodes_.begin() + p_start, p.nodes_.end());
}

void Path::AppendPathWithGap(const Path& p, int gap_length, int p_start) {
  nodes_.push_back(MakeGap(gap_length));
  nodes_.insert(nodes_.end(), p.nodes_.begin() + p_start, p.nodes_.end());
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

string Path::ToString(bool with_endings) const {
  string ret = "";
  if (with_endings) {
    assert((int) nodes_[0]->str_.size() > nodes_[0]->graph_->k_ - 1);
    ret = ReverseSeq(nodes_[0]->rc_->str_).substr(0, nodes_[0]->graph_->k_ - 1);
  }
  for (auto &n: nodes_) {    
    ret += n->str_;
  }
  return ret;
}

bool Path::IsSame(const Path& p) const {
  if (p.nodes_.size() != nodes_.size()) return false;
  if (p.nodes_ == nodes_) return true;
  for (size_t i = 0; i < nodes_.size(); i++) {
    if (nodes_[i]->rc_ != p.nodes_[nodes_.size() - 1 - i]) return false;
  }
  return true;
}

bool Path::ExtendRandomly(int big_node_threshold, int step_threshold, int distance_threshold) {
  int added_distance = 0;
  int added_steps = 0;
  do {
    Node* last_node = nodes_.back();
    if (last_node->next_.size() == 0) return false;
    Node* next_node = last_node->next_[rand()%last_node->next_.size()];
    nodes_.push_back(next_node);
    if ((int)next_node->str_.size() >= big_node_threshold) {
      return true;
    }
    added_steps += 1;
    added_distance += next_node->str_.size();
  } while (added_distance <= distance_threshold && added_steps <= step_threshold);
  return false;
}

void PathsToFasta(const vector<Path>& paths, ostream &of) {
  for (auto &p: paths) {
    of << ">" << p.ToDebugString() << endl;
    of << p.ToString(true) << endl;
  }
}

void ComparePathSets(const vector<Path>& a,
                     const vector<Path>& b,
                     vector<Path>& added,
                     vector<Path>& removed) {
  for (auto &pb: b) {
    bool found = false;
    for (auto &pa: a) {
      if (pa.IsSame(pb)) {
        found = true;
        break;
      }
    }
    if (!found) {
      added.push_back(pb);
    }
  }
  for (auto &pa: a) {
    bool found = false;
    for (auto &pb: b) {
      if (pa.IsSame(pb)) {
        found = true;
        break;
      }
    }
    if (!found) {
      removed.push_back(pa);
    }
  }
}
