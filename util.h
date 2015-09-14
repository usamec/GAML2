#ifndef UTIL_H__
#define UTIL_H__

#include <string>

using namespace std;

inline char ReverseBase(char c) {
  if (c == 'A') return 'T';
  if (c == 'C') return 'G';
  if (c == 'G') return 'C';
  if (c == 'T') return 'A';
  return 'N';
}

inline string ReverseSeq(const string& s) {
  string ret;
  for (int i = s.size() - 1; i >= 0; i--) {
    ret += ReverseBase(s[i]);
  }
  return ret;
}

#endif
