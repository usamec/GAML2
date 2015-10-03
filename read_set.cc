#include "read_set.h"
#include "hash_util.h"
#include "util.h"
#include <algorithm>
#include <deque>
#include <unordered_set>

void StandardReadIndex::AddRead(int id, const string& data) {
  for (size_t i = 0; i + k_ <= data.size(); i++) {
    index_[data.substr(i, k_)].push_back(make_pair(id, i));
  }
}

vector<CandidateReadPosition> StandardReadIndex::GetReadCandidates(const string& genome) const {
  //read_id, diagonal / 5
  static unordered_set<pair<int, int>> found_cands;
  found_cands.clear();
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

template<class TIndex>
vector<ReadAlignment> ReadSet<TIndex>::GetAlignments(const string& genome) const {
  vector<ReadAlignment> ret;
  GetAlignments(genome, false, ret);
  string reversed_genome = ReverseSeq(genome);
  GetAlignments(reversed_genome, true, ret);
  return ret;
}

template<class TIndex>
void ReadSet<TIndex>::GetAlignments(const string& genome,
                                    bool reversed,
                                    vector<ReadAlignment>& output) const {
  vector<CandidateReadPosition> candidates = index_.GetReadCandidates(genome);

  sort(candidates.begin(), candidates.end());

  int last_read_id = -1;
  vector<ReadAlignment> buffer;
  for (auto &cand: candidates) {
    if (cand.read_id != last_read_id) {
      output.insert(output.end(), buffer.begin(), buffer.end());
      buffer.clear();
    }
    last_read_id = cand.read_id;
    ReadAlignment al;
    if (ExtendAlignment(cand, genome, al)) {
      if (reversed) {
        al.genome_pos = genome.size() - al.genome_pos - reads_[cand.read_id].size();
      }
      auto it = find_if(
          buffer.begin(), buffer.end(),
          [&al](const ReadAlignment& a) { return a.read_id == al.read_id && 
                                                 a.genome_pos == al.genome_pos; });
      if (it == buffer.end()) {
        al.reversed = reversed;
        buffer.push_back(al);
      } else {
        it->dist = min(it->dist, al.dist);
      }
    }
  }
  output.insert(output.end(), buffer.begin(), buffer.end());
}

template<class TIndex>
void ReadSet<TIndex>::VisitedPositions::Prepare(int max_err, int read_size) {
  vp_.resize(max_err + 2);
  for (auto &x: vp_) {
    x.resize(read_size+2);
    for (auto &y: x) {
      y.clear();
    }
  }
}


template<class TIndex>
bool ReadSet<TIndex>::VisitedPositions::IsVisited(
    pair<int, pair<int, int>> pos) const {
  pos.second.first += 1;
  return find(vp_[pos.first][pos.second.first].begin(),
              vp_[pos.first][pos.second.first].end(),
              pos.second.second) != vp_[pos.first][pos.second.first].end();
}

template<class TIndex>
void ReadSet<TIndex>::VisitedPositions::Add(pair<int, pair<int, int>> pos) {
  pos.second.first += 1;
  vp_[pos.first][pos.second.first].push_back(pos.second.second);
}


template<class TIndex>
bool ReadSet<TIndex>::ExtendAlignment(const CandidateReadPosition& candidate,
                                      const string& genome,
                                      ReadAlignment& al) const {
  int max_err_start = 6;
  int max_err = max_err_start;
  // Static - we reuse memory and make fewer allocations
  static VisitedPositions visited_positions;

  // indexing: distance -> read_pos -> list of genome_positions
  auto& read = reads_[candidate.read_id];
  int total_errs = 0;

  visited_positions.Prepare(max_err, read.size());

  // distance, (read_pos, genome_pos)
  static deque<pair<int, pair<int, int>>> fr;
  fr.clear();
  fr.push_back(make_pair(0, make_pair(candidate.read_pos+1, candidate.genome_pos+1)));

  while (!fr.empty()) {
    auto x = fr.front();
    fr.pop_front();

    if (visited_positions.IsVisited(x)) continue;
    visited_positions.Add(x);

    if (x.second.first == read.size()) {
      total_errs = x.first;
      break;
    }

    if (x.first > max_err) {
      return false; 
    }

    if (x.second.second < genome.size()) {
      // match / mismatch
      if (genome[x.second.second] == read[x.second.first]) {
        fr.push_front(make_pair(x.first, make_pair(x.second.first+1, x.second.second+1)));
        // if matches we greedily continue
        continue;
      }
      // mismatch
      fr.push_back(make_pair(x.first+1, make_pair(x.second.first+1, x.second.second+1)));

      // delete from genome
      fr.push_back(make_pair(x.first+1, make_pair(x.second.first, x.second.second+1)));
    }
    // insert to genome
    fr.push_back(make_pair(x.first+1, make_pair(x.second.first+1, x.second.second)));
  }

  max_err -= total_errs;

  fr.clear();
  fr.push_back(make_pair(0, make_pair(candidate.read_pos, candidate.genome_pos)));
  while (!fr.empty()) {
    auto x = fr.front();
    fr.pop_front();

    if (visited_positions.IsVisited(x)) continue;
    visited_positions.Add(x);

    if (x.second.first == -1) {
      al.read_id = candidate.read_id; 
      al.dist = total_errs + x.first;
      al.genome_pos = x.second.second + 1;
      return true;
    }

    if (x.first > max_err) {
      return false;
    }

    if (x.second.second >= 0) {
      // match / mismatch
      if (genome[x.second.second] == read[x.second.first]) {
        fr.push_front(make_pair(x.first, make_pair(x.second.first-1, x.second.second-1)));
        // if matches we greedily continue
        continue;
      }
      // mismatch
      fr.push_back(make_pair(x.first+1, make_pair(x.second.first-1, x.second.second-1)));

      // delete from genome
      fr.push_back(make_pair(x.first+1, make_pair(x.second.first, x.second.second-1)));
    }
    // insert to genome
    fr.push_back(make_pair(x.first+1, make_pair(x.second.first-1, x.second.second)));
  }
  return false;
}

template class ReadSet<StandardReadIndex>;
