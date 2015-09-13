#ifndef HASH_UTIL_H__
#define HASH_UTIL_H__

#include <unordered_set>

using namespace std;

template <class T>
inline void hash_combine(std::size_t & seed, const T & v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
  template<typename S, typename T> struct hash<pair<S, T>> {
    inline size_t operator()(const pair<S, T> & v) const {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };

  template<typename T> struct hash<vector<T>> {
    inline size_t operator()(const vector<T>& v) const {
      size_t seed = 0;
      for (int i = 0; i < v.size(); i++) {
        ::hash_combine(seed, v[i]);
      }
      return seed;
    }
  };
}

#endif
