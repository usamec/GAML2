#include "read_set.h"
#include "hash_util.h"
#include <unordered_set>

void StandardReadIndex::AddRead(int id, const string& data) {
  for (size_t i = 0; i + k_ <= data.size(); i++) {
    index_[data.substr(i, k_)].push_back(make_pair(id, i));
  }
}

vector<CandidateReadPosition> StandardReadIndex::GetReadCandidates(const string& genome) const {
  //read_id, diagonal / 5
  unordered_set<pair<int, int>> found_cands;
  vector<CandidateReadPosition> ret;

  for (size_t i = 0; i + k_ <= genome.size(); i++) {
    auto it = index_.find(genome.substr(i, k_));
    if (it == index_.end()) continue;
    for (auto &e: it->second) {
      int coord = (i - e.second) / 5;
      if (found_cands.count(make_pair(e.first, coord))) {
        continue;
      }
      found_cands.insert(make_pair(e.first, coord));
      ret.push_back(CandidateReadPosition(e.first, i, e.second));
    }
  }
  return ret;
}

template<class TIndex>
void ReadSet<TIndex>::LoadReadSet(istream& is) {
  string l1, l2, l3, l4;
  int id = 0;
  while (getline(is, l1)) {
    getline(is, l2);
    getline(is, l3);
    getline(is, l4);
    reads_.push_back(l2);
    index_.AddRead(id, l2);
    id++;
  }
}

template class ReadSet<StandardReadIndex>;
