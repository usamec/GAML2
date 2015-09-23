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

  ASSERT_EQ(1, rs.size());
  EXPECT_EQ("AAGGCTG", rs[0]);
}

TEST(ReadSetTest, ExtendAlignTest) {
  stringstream ss;
  ss << "@a" << endl;
  ss << "AAAACCCCTTTTGGGG" << endl;
  ss << "+" << endl;
  ss << "AAAAAAAAAAAAAAAA" << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss);

  CandidateReadPosition candidate(0, 0, 0);
  ReadAlignment al;
  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "AAAACCCCTTTTGGGGAAA", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(0, al.dist);
  EXPECT_EQ(0, al.genome_pos);
  
  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "AAAACCCCTTTTGGGG", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(0, al.dist);
  EXPECT_EQ(0, al.genome_pos);
  
  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "AAAACCCCTTTTGGTG", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(1, al.dist);
  EXPECT_EQ(0, al.genome_pos);

  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "AAAACCCCTTTTGGT", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(2, al.dist);
  EXPECT_EQ(0, al.genome_pos);

  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "AAAACCAAAATTGGGG",  al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(4, al.dist);
  EXPECT_EQ(0, al.genome_pos);

  EXPECT_EQ(false, rs.ExtendAlignment(candidate, "AAAACCAAAAAACGGG",  al));
  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "AAAACCAAAAAAGGGG",  al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(6, al.dist);
  EXPECT_EQ(0, al.genome_pos);


  candidate.genome_pos = 3;
  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "CACAAAACCCCTTTTGGGG", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(0, al.dist);
  EXPECT_EQ(3, al.genome_pos);
  candidate.genome_pos = 5;
  candidate.read_pos = 2;
  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "CACAAAACCCCTTTTGGGG", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(0, al.dist);
  EXPECT_EQ(3, al.genome_pos);

  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "CACTTAACCCCTTTTGGGG", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(2, al.dist);
  EXPECT_EQ(3, al.genome_pos);

  EXPECT_EQ(true, rs.ExtendAlignment(candidate, "CACATAACCCCTTTTGGGG", al));
  EXPECT_EQ(0, al.read_id);
  EXPECT_EQ(1, al.dist);
  EXPECT_EQ(3, al.genome_pos);
}
