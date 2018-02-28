#include "graph.h"
#include <boost/algorithm/string.hpp>
#include <queue>
#include <unordered_set>
#include <iostream>

using boost::is_any_of;
using boost::split;

int VelvetNumToMyNum(int x) {
  if (x > 0) {
    return (x-1)*2;
  } else {
    return (-x-1)*2+1;
  }
}

pair<int, int> LoadArc(istream& is) {
  string arc_str;
  if (!getline(is, arc_str)) {
    return make_pair(-1, -1);
  }
  if (arc_str.substr(0, 3) != "ARC") {
    return make_pair(-1, -1);
  }
 
  vector<string> arc_tokens;
  split(arc_tokens, arc_str, is_any_of("\t"));

  return make_pair(VelvetNumToMyNum(stoi(arc_tokens[1])),
                   VelvetNumToMyNum(stoi(arc_tokens[2])));
}

Graph* LoadGraph(const string &filename) {
  ifstream is(filename);
  if (!is) return NULL;
  return LoadGraph(is);
}

Graph* LoadGraph(istream &is) {
  string header;
  getline(is, header);
  vector<string> header_parts;
  split(header_parts, header, is_any_of("\t"));
  Graph* g = new Graph();
  int n_nodes = stoi(header_parts[0]);
  g->k_ = stoi(header_parts[2]);
  for (int i = 0; i < n_nodes; i++) {
    Node *a, *b;
    tie(a, b) = LoadNode(is, i, g);
    g->nodes_.push_back(a);
    g->nodes_.push_back(b);
  }

  while (is) {
    int a, b;
    tie(a, b) = LoadArc(is);
    if (a != -1) {
      g->nodes_[a]->AddNext(g->nodes_[b]);
      g->nodes_[b]->rc_->AddNext(g->nodes_[a]->rc_);
    }
  }
  return g;
}

vector<Node*> Graph::GetBigNodes(int threshold) const {
  vector<Node*> ret;
  for (size_t i = 0; i < nodes_.size(); i+= 2) {
    if (nodes_[i]->IsBig(threshold)) {
      ret.push_back(nodes_[i]);
    }
  }
  return ret;
}


// BFS but ignoring big nodes (based on threshold)?
// + ignoring reversed complement nodes
vector<Node*> Graph::ReachForwardWithThreshold(Node* start, int threshold) const {
  unordered_set<int> visited;
  queue<Node*> fr;
  fr.push(start);
  visited.insert(start->id_);
  vector<Node*> ret;

  while (!fr.empty()) {
    Node *x = fr.front();
    fr.pop();
    ret.push_back(x);
    for (auto &nx: x->next_) {
      if (nx->IsBig(threshold)) {
        continue;
      }
      // unordered_set::count := # of occurences in the set (may be 0 or 1)
      if (visited.count(nx->id_)) {
        continue;
      }
      visited.insert(nx->id_);
      fr.push(nx);
    }
  }

  return ret;
}

// BFS but ignoring big nodes (based on threshold)?
vector<Node*> Graph::ReachLocalWithThreshold(Node* start, int threshold) const {
  unordered_set<int> visited;
  queue<Node*> fr;
  fr.push(start);
  visited.insert(start->id_);
  vector<Node*> ret;

  while (!fr.empty()) {
    Node *x = fr.front();
    fr.pop();
    ret.push_back(x);
    for (auto &nx: x->next_) {
      if (nx->IsBig(threshold)) {
        continue;
      }
      if (visited.count(nx->id_)) {
        continue;
      }
      visited.insert(nx->id_);
      fr.push(nx);
    }
    for (auto &nxr: x->rc_->next_) {
      auto &nx = nxr->rc_;
      if (nx->IsBig(threshold)) {
        continue;
      }
      if (visited.count(nx->id_)) {
        continue;
      }
      visited.insert(nx->id_);
      fr.push(nx);

    }
  }

  return ret;
}

// https://www.wikiwand.com/en/Drainage_basin
// BFS on in-endges
// get nodes for with the given node is reachable
vector<Node*> Graph::GetDrainageBasinForNode(Node* target) {
  unordered_set<int> visited_ids;
  queue<Node*> Q;
  Q.push(target);
  visited_ids.insert(target->id_);

  vector<Node*> ret;

  while (!Q.empty()) {
    Node* x = Q.front();
    Q.pop();
    ret.push_back(x);

    for (auto &px: x->prev_) {
      if (visited_ids.count(px->id_) > 0) continue;
      Q.push(px);
      visited_ids.insert(px->id_);
      }
    }
  return ret;
}