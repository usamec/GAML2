#include "path.h"
#include "graph.h"
#include "util.h"
#include <cassert>
#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <iostream>

vector<Path> BuildPathsFromSingleNodes(const vector<Node*>& nodes) {
  vector<Path> ret;
  for (auto &n: nodes) {
    ret.push_back(Path({n}, PATH_BIGNODE));
  }
  return ret;
}

// nodes have to be adjacent (in the same direction) in the graph
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

string PathsToDebugStringShort(const vector<Path>& paths) {
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

string PathsToDebugString(const vector<Path>& paths) {
  stringstream ret;
  ret << paths.size() << " paths:\n";
  for (size_t i = 0; i < paths.size(); i++) {
    ret << "[" << paths[i].history_ <<  "] " << paths[i].ToDebugString();
    if (i + 1 != paths.size()) {
      ret << "\n";
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

Path Path::GetReverse() const {
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

bool Path::IsSameNoReverse(const Path& p) const {
  if (p.nodes_.size() != nodes_.size()) return false;
  if (p.nodes_ == nodes_) return true;
  return false;
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

Path Path::CutAt(int pos, int big_node_threshold) {
  int part1_end = pos - 1;
  while (part1_end >= 0 && !nodes_[part1_end]->IsBig(big_node_threshold)) part1_end--;
  int part2_start = pos;
  while (part2_start < nodes_.size() && !nodes_[part2_start]->IsBig(big_node_threshold)) part2_start++;
  Path p2;
  if (part2_start < nodes_.size()) p2 = Path(vector<Node*>(nodes_.begin() + part2_start, nodes_.end()), history_);

  nodes_ = vector<Node*>(nodes_.begin(), nodes_.begin() + (part1_end + 1));
  return p2;
}
bool Path::isDisjoint(const Path &p) const {
  // complementary nodes also not allowed
  unordered_set<int> nodes_ids;
  for (auto node: nodes_) {
    nodes_ids.insert(node->id_);
  }
  for (auto node: p.nodes_) {
    if (nodes_ids.count(node->id_) > 0) return false;
    // check for complementar
    if (nodes_ids.count(node->rc_->id_) > 0) return false;
  }
  return true;
}
bool Path::isPartlyDisjoint(const Path &p, int dist) const {
  // complementary nodes also not allowed
  unordered_set<int> nodes_ids;

  int len = 0;
  for (size_t i = 0; i < nodes_.size() && len < dist; i++) {
    const auto &node = nodes_[i];
    nodes_ids.insert(node->id_);
    len += node->str_.size();
  }
  len = 0;
  for (size_t i = nodes_.size() - 1; i >= 0 && len < dist; i--) {
    const auto &node = nodes_[i];
    nodes_ids.insert(node->id_);
    len += node->str_.size();
  }

  // checking begin and end of path p
  len = 0;
  for (size_t i = 0; i < p.nodes_.size() && len < dist; i++) {
    const auto &node = p.nodes_[i];
    if (nodes_ids.count(node->id_) > 0 || nodes_ids.count(node->rc_->id_) > 0) return false;
    len += node->str_.size();
  }
  len = 0;
  for (size_t i = p.nodes_.size() - 1; i >= 0 && len < dist; i--) {
    const auto &node = p.nodes_[i];
    if (nodes_ids.count(node->id_) > 0 || nodes_ids.count(node->rc_->id_) > 0) return false;
    len += node->str_.size();
  }

  return true;
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
