#include "read_set.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>

TEST(ReadSetTest, LoadTest) {
  stringstream ss;
  ss << "@aaa" << endl;
  ss << "AAGGCTG" << endl;
  ss << "+" << endl;
  ss << "AAAAAAA" << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss);
}
