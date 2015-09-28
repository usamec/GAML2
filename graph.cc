#include "graph.h"
#include <boost/algorithm/string.hpp>

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
