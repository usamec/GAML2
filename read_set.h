#ifndef READ_SET_H__
#define READ_SET_H__

#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <gtest/gtest.h>
#include "Sequence.h"
#include "DalignWrapper.h"
#include <unordered_set>
using namespace std;

struct CandidateReadPosition {
  CandidateReadPosition() {};
  CandidateReadPosition(int read_id_, int genome_pos_, int read_pos_) :
      read_id(read_id_), genome_pos(genome_pos_), read_pos(read_pos_) {}

  int read_id, genome_pos, read_pos;
};

inline bool operator==(const CandidateReadPosition& a, const CandidateReadPosition& b) {
  return a.read_id == b.read_id && a.genome_pos == b.genome_pos && a.read_pos == b.read_pos;
}

inline bool operator<(const CandidateReadPosition& a, const CandidateReadPosition& b) {
  if (a.read_id < b.read_id) return true;
  if (a.read_id > b.read_id) return false;
  if (a.genome_pos < b.genome_pos) return true;
  if (a.genome_pos > b.genome_pos) return false;
  return a.read_pos < b.read_pos;
}

struct ReadAlignment {
  ReadAlignment() {}
  ReadAlignment(int read_id_, int genome_pos_, int dist_, bool reversed_) :
      read_id(read_id_), genome_pos(genome_pos_), dist(dist_), reversed(reversed_) {}

  int read_id;
  int genome_pos;
  int dist;
  bool reversed;
};

inline bool operator<(const ReadAlignment& a, const ReadAlignment& b) {
  return a.read_id < b.read_id;
}

class StandardReadIndex {
 public:
  StandardReadIndex(int k = 13): k_(k) {}
  void AddRead(int id, const string& data);

  vector<CandidateReadPosition> GetReadCandidates(const string& genome) const;

  int k_;
  // (read_id, pos_in_read)
  unordered_map<string, vector<pair<int,int>>> index_;
};

class RandomIndex {
 public:
  RandomIndex(int k = 13): k_(k) {}
  void AddRead(int id, const string& data);

  vector<CandidateReadPosition> GetReadCandidates(const string& genome) const;

  int k_;
  // (read_id, pos_in_read)
  unordered_map<string, vector<pair<int,int>>> index_; 
};

template<class TIndex=RandomIndex>
class ReadSet {
  class VisitedPositions {
   vector<vector<int>> vp_;
   int offset_;
   int offset2_;
   int cur_iter_;
   public:
    VisitedPositions() : cur_iter_(0) {}
    void Prepare(int offset, int read_size);
    
    bool IsVisited(pair<int, pair<int, int>> pos) const {
      return vp_[pos.second.first + 10][pos.second.second - offset_ + offset2_ + 10] == cur_iter_;
    }
    void Add(pair<int, pair<int, int>> pos) {
      vp_[pos.second.first + 10][pos.second.second - offset_ + offset2_ + 10] = cur_iter_;
    }
  };

 public:
  ReadSet() {}

  void LoadReadSet(const string& filename) {
    ifstream is(filename);
    LoadReadSet(is);
  }

  void LoadReadSet(istream& is);

  // Two sided get
  vector<ReadAlignment> GetAlignments(const string& genome) const;

  size_t size() const {
    return reads_.size();
  }

  const string& operator[](int i) const {
    return reads_[i];
  }

 private:
  // One sided get
  void GetAlignments(const string& genome, bool reversed, vector<ReadAlignment>& output) const;

  bool ExtendAlignment(const CandidateReadPosition& candidate, const string& genome,
                       ReadAlignment& al) const;

  vector<string> reads_;
  TIndex index_;

  FRIEND_TEST(ReadSetTest, ExtendAlignTest);
};


struct ReadAlignmentPacBio {
  ReadAlignmentPacBio() {}
  ReadAlignmentPacBio(int read_id_, int genome_first_, int genome_last_, int read_first_, int read_last_, int dist_, bool reversed_) :
      read_id(read_id_), genome_first(genome_first_), genome_last(genome_last_), read_first(read_first_), 
      read_last(read_last_), dist(dist_), reversed(reversed_) {}

  int read_id;
  int genome_first;
  int genome_last;
  int read_first;
  int read_last;
  int dist;
  bool reversed;
};

template<class TIndex=RandomIndex>
class ReadSetPacBio {
  
  class AlignedPairsSet {
    unordered_set<long long> alignedPairsSet;
  public:
      void MarkAsAligned(vector<pair<int, int>> &alignedPairsVector);
      bool IsAligned(pair<int, int> seed);
  };
  
public:
  ReadSetPacBio() {
    SetParameters(0.7, {0.25, 0.25, 0.25, 0.25}, 100);
  }
    
  void LoadReadSet(const string& filename);

  void LoadReadSet(istream& is);

  // Two sided get
  vector<ReadAlignmentPacBio> GetAlignments(const string& genome);

  size_t size() const {
    return reads_.size();
  }

  const string& operator[](int i) const {
    return reads_[i].GetData();
  }
  
  void SetParameters(float corelation, const array<float, 4>& frequencies, int minSufficientLength_);
  
private:
  
  // One sided get
  void GetAlignments(Sequence& genome, bool reversed, vector<ReadAlignmentPacBio>& output);
  
  vector<Sequence> reads_;
  
  TIndex index_;
  
  DalignWrapper dalign_;
  
  int minSufficientLength;
};

#endif
