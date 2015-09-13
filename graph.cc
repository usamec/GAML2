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
  getline(is, arc_str);
 
  vector<string> arc_tokens;
  split(arc_tokens, arc_str, is_any_of("\t"));

  return make_pair(VelvetNumToMyNum(stoi(arc_tokens[1])),
                   VelvetNumToMyNum(stoi(arc_tokens[2])));
}
