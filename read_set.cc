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

void RandomIndex::AddRead(int id, const string& data) {
  if ((int) data.size() < k_) return;
  for (int i = 0; i < 3; i++) {
    int p = rand()%(data.size() - k_ + 1);
    index_[data.substr(p, k_)].push_back(make_pair(id, p));
  }
}

vector<CandidateReadPosition> RandomIndex::GetReadCandidates(const string& genome) const {
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
    if (id % 10000 == 0) {
      printf("\rLoaded %d reads", id);
      fflush(stdout);
    }
  }
  printf("\n");
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
void ReadSet<TIndex>::VisitedPositions::Prepare(int offset, int read_size) {
  if (read_size * 2 + 40 > (int)vp_.size()) {
    vp_.resize(read_size * 2 + 40);
    for (auto &x: vp_) {
      x.resize(read_size * 2 + 40);
    }
  }
  offset_ = offset;
  offset2_ = read_size;
  cur_iter_++;
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

  visited_positions.Prepare(candidate.genome_pos, read.size());

  // distance, (read_pos, genome_pos)
  static deque<pair<int, pair<int, int>>> fr;
  fr.clear();
  fr.push_back(make_pair(0, make_pair(candidate.read_pos+1, candidate.genome_pos+1)));

  while (!fr.empty()) {
    auto x = fr.front();
    fr.pop_front();

    if (x.second.first == (int)read.size()) {
      total_errs = x.first;
      break;
    }

    if (x.first > max_err) {
      return false; 
    }

    if (x.second.second < (int)genome.size()) {
      // match / mismatch
      if (genome[x.second.second] == read[x.second.first]) {
        auto nx = make_pair(x.first, make_pair(x.second.first+1, x.second.second+1));
        fr.push_front(nx);
        visited_positions.Add(nx);
        // if matches we greedily continue
        continue;
      }
      // mismatch
      {
        auto nx = make_pair(x.first+1, make_pair(x.second.first+1, x.second.second+1));
        if (!visited_positions.IsVisited(nx)) {
          fr.push_back(nx);
          visited_positions.Add(nx);
        }
      }

      // delete from genome
      {
        auto nx = make_pair(x.first+1, make_pair(x.second.first, x.second.second+1));
        if (!visited_positions.IsVisited(nx)) {
          fr.push_back(nx);
          visited_positions.Add(nx);
        }
      }
    }
    // insert to genome
    {
      auto nx = make_pair(x.first+1, make_pair(x.second.first+1, x.second.second));
      if (!visited_positions.IsVisited(nx)) {
        fr.push_back(nx);
        visited_positions.Add(nx);
      }
    }
  }

  max_err -= total_errs;

  fr.clear();
  fr.push_back(make_pair(0, make_pair(candidate.read_pos, candidate.genome_pos)));
  while (!fr.empty()) {
    auto x = fr.front();
    fr.pop_front();

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
        auto nx = make_pair(x.first, make_pair(x.second.first-1, x.second.second-1));
        fr.push_front(nx);
        visited_positions.Add(nx);
        // if matches we greedily continue
        continue;
      }
      // mismatch
      {
        auto nx = make_pair(x.first+1, make_pair(x.second.first-1, x.second.second-1));
        if (!visited_positions.IsVisited(nx)) {
          fr.push_back(nx);
          visited_positions.Add(nx);
        }
      }

      // delete from genome
      {
        auto nx = make_pair(x.first+1, make_pair(x.second.first, x.second.second-1));
        if (!visited_positions.IsVisited(nx)) {
          fr.push_back(nx);
          visited_positions.Add(nx);
        }
      }
    }
    // insert to genome
    {
      auto nx = make_pair(x.first+1, make_pair(x.second.first-1, x.second.second));
      if (!visited_positions.IsVisited(nx)) {
        fr.push_back(nx);
        visited_positions.Add(nx);
      }
    }
  }
  return false;
}


template<class TIndex>
void ReadSetPacBio<TIndex>::AlignedPairsSet::MarkAsAligned(vector<pair<int, int> >& alignedPairsVector) {
  long long previous = -1;
  for (auto &p : alignedPairsVector) {
    long long numValue = (((long long)p.first >> 3) << 32) + (p.second >> 3);
    if (numValue != previous) {
      alignedPairsSet.insert(numValue);
      previous = numValue;
    }
  }
}


template<class TIndex>
bool ReadSetPacBio<TIndex>::AlignedPairsSet::IsAligned(pair<int,int> seed) {
  long long numValue = (((long long)(seed.first) >> 3) << 32) + ((seed.second) >> 3);
  return alignedPairsSet.find(numValue) != alignedPairsSet.end();
}


template<class TIndex>
void ReadSetPacBio<TIndex>::LoadReadSet(istream& is) {
  string l1, l2, l3, l4;
  int id = 0;
  while (getline(is, l1)) {
    getline(is, l2);
    getline(is, l3);
    getline(is, l4);
    reads_.push_back(Sequence(l2));
    index_.AddRead(id, l2);
    id++;
    if (id % 10000 == 0) {
      printf("\rLoaded %d reads", id);
      fflush(stdout);
    }
  }
  printf("\n");
}

template<class TIndex>
void ReadSetPacBio<TIndex>::SetParameters(float corelation, const array<float, 4>& frequencies, int minSufficientLength_) {
  dalign_.SetAligningParameters(corelation, 50, frequencies);
  minSufficientLength = minSufficientLength_;
}


template<class TIndex>
vector<ReadAlignmentPacBio> ReadSetPacBio<TIndex>::GetAlignments(const string& genome) {
  vector<ReadAlignmentPacBio> result;
  Sequence forwardGenome(genome);
  Sequence backwardGenome(forwardGenome, true);
  
  GetAlignments(forwardGenome, false, result);
  GetAlignments(backwardGenome, true, result);
  
  return result;
}

template<class TIndex>
void ReadSetPacBio<TIndex>::GetAlignments(Sequence& genome, bool reversed, vector<ReadAlignmentPacBio>& output) {
  vector<CandidateReadPosition> candidates = index_.GetReadCandidates(genome.GetData());
  
  sort(candidates.begin(), candidates.end());
  AlignedPairsSet alignedPairs;
  
  int lastId = -1;
  for (auto& candidate : candidates) {
    if (candidate.read_id != lastId) {
      alignedPairs = AlignedPairsSet();
    }
    lastId = candidate.read_id;
    
    if (alignedPairs.IsAligned(make_pair(candidate.genome_pos, candidate.read_pos))) {
      continue;
    }
    
    Alignment al;
    Sequence &read = reads_[candidate.read_id];
    dalign_.ComputeAlignment(genome, read, pair<int, int>(candidate.genome_pos, candidate.read_pos), al);
    int length = (al.GetLengthOnA() + al.GetLengthOnB()) / 2;
    
    if (length >= minSufficientLength) {
      
      al.ComputeTrace();
      
      vector<pair<int, int>> pairs;
      al.GetAlignedPairs(pairs);
      alignedPairs.MarkAsAligned(pairs);
      
      ReadAlignmentPacBio newAlignment;
      newAlignment.read_id = candidate.read_id;
      newAlignment.read_first = al.GetPosOnB().first;
      newAlignment.read_last = al.GetPosOnB().second;
      newAlignment.reversed = reversed;
      newAlignment.dist = al.GetNumDifferences();
      if (reversed) {
        newAlignment.genome_first = genome.GetData().length() - al.GetPosOnA().second;
        newAlignment.genome_last = genome.GetData().length() - al.GetPosOnA().first;
      } else {
        newAlignment.genome_first = al.GetPosOnA().first;
        newAlignment.genome_last = al.GetPosOnA().second;
      }
      output.push_back(newAlignment);
    }
  }
}

template class ReadSet<StandardReadIndex>;
template class ReadSet<RandomIndex>;
template class ReadSetPacBio<StandardReadIndex>;
template class ReadSetPacBio<RandomIndex>;
