#ifndef READ_SET_H__
#define READ_SET_H__

#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>

using namespace std;

struct CandidateReadPosition {
  CandidateReadPosition() {};
  CandidateReadPosition(int read_id_, int genome_pos_, int read_pos_) :
      read_id(read_id_), genome_pos(genome_pos_), read_pos(read_pos_) {}

  int read_id, genome_pos, read_pos;
};

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
 public:
  ReadSet() {}

  void LoadReadSet(const string& filename) {
    ifstream is(filename);
    LoadReadSet(is);
  }

  void LoadReadSet(istream& is);

  vector<string> reads_;
  TIndex index_;
};

#endif
