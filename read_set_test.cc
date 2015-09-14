#include "read_set.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>

TEST(StandardReadIndexTest, GetCandidatesTest) {
  StandardReadIndex index(13);
  //             12345678901234567890
  index.AddRead(1, "AAAAAAAAAAAAAACCCTTT");
  index.AddRead(2, "AAAAAACCCTTTTTTTTTTT");

  vector<CandidateReadPosition> result = index.GetReadCandidates(
      "AAAAACCCTTTTTTTTTTTTTTTTTTTTTTTTT");
  vector<CandidateReadPosition> expected;
  expected.push_back(CandidateReadPosition(2, 0, 1));
  EXPECT_EQ(expected, result);

  expected.clear();
  EXPECT_EQ(expected, index.GetReadCandidates("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"));
}

TEST(ReadSetTest, LoadTest) {
  stringstream ss;
  ss << "@aaa" << endl;
  ss << "AAGGCTG" << endl;
  ss << "+" << endl;
  ss << "AAAAAAA" << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss);
}
