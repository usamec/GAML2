#ifndef READ_SET_H__
#define READ_SET_H__

#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <gtest/gtest.h>

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

class StandardReadIndex {
 public:
  StandardReadIndex(int k = 13): k_(k) {}
  void AddRead(int id, const string& data);

  vector<CandidateReadPosition> GetReadCandidates(const string& genome) const;

  int k_;
  // (read_id, pos_in_read)
  unordered_map<string, vector<pair<int,int>>> index_;
};

template<class TIndex=StandardReadIndex>
class ReadSet {
  class VisitedPositions {
   vector<vector<vector<int>>> vp_;
   public:
    void Prepare(int max_err, int read_size);
    
    bool IsVisited(pair<int, pair<int, int>> pos) const;
    void Add(pair<int, pair<int, int>> pos);
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

#endif
