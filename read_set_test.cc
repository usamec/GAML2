#include "read_set.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>
#include "util.h"

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

TEST(ReadSetTest, GetAlignmentsTest) {
  stringstream ss;
  ss << "@a" << endl;
  ss << "AAAACCCCTTTTGGGG" << endl;
  ss << "+" << endl;
  ss << "AAAAAAAAAAAAAAAA" << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss);

  vector<ReadAlignment> als = rs.GetAlignments("AAAACCCCTTTTGGGG");
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);
  EXPECT_EQ(false, als[0].reversed);
  EXPECT_EQ(0, als[0].genome_pos);

  als = rs.GetAlignments("AAAACCCCTTTTGGCG");
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(1, als[0].dist);
  EXPECT_EQ(false, als[0].reversed);
  EXPECT_EQ(0, als[0].genome_pos);

  als = rs.GetAlignments("CCCCAAAACCCCTTTTGGCG");
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(1, als[0].dist);
  EXPECT_EQ(false, als[0].reversed);
  EXPECT_EQ(4, als[0].genome_pos);

  als = rs.GetAlignments("CCCCAAAAGGGGTTTT");
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);
  EXPECT_EQ(true, als[0].reversed);
  EXPECT_EQ(0, als[0].genome_pos);

  als = rs.GetAlignments("CCCCAAAAGGGGTTTTAA");
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);
  EXPECT_EQ(true, als[0].reversed);
  EXPECT_EQ(0, als[0].genome_pos);

  als = rs.GetAlignments("AAACCCCAAAAGGGGTTTT");
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);
  EXPECT_EQ(true, als[0].reversed);
  EXPECT_EQ(3, als[0].genome_pos);

  als = rs.GetAlignments("AAACCCCAAAAGGGGTTTTGG");
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);
  EXPECT_EQ(true, als[0].reversed);
  EXPECT_EQ(3, als[0].genome_pos);
}

TEST(ReadSetTest, GetAlignmentsTestBig) {
  string part1 = "";
  string part2 = "";
  srand(47);
  char alph[] = "ACGT";
  for (int i = 0; i < 50; i++) {
    part1 += alph[rand()%4];
    part2 += alph[rand()%4];
  }

  stringstream ss;
  ss << "@a" << endl;
  ss << part1 << part2 << endl;
  ss << "+" << endl;
  ss << part1 << part2 << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss);

  string genome = "ACGTTT" + part1 + part2 + "GTCT";
  vector<ReadAlignment> als = rs.GetAlignments(genome);
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(false, als[0].reversed);
  EXPECT_EQ(6, als[0].genome_pos);
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);

  genome = "ACGTTT" + ReverseSeq(part1+part2) + "GTCT";
  als = rs.GetAlignments(genome);
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(true, als[0].reversed);
  EXPECT_EQ(6, als[0].genome_pos);
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);

  genome = "ACGTTT" + part1 + "ACTGAA" + part2 + "GTCT";
  als = rs.GetAlignments(genome);
  ASSERT_EQ(1, als.size());
  EXPECT_EQ(false, als[0].reversed);
  EXPECT_EQ(6, als[0].genome_pos);
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(6, als[0].dist);
}

TEST(ReadSetTest, GetAlignmentsTestMoreReads) {
  string part1 = "";
  string part2 = "";
  string part3 = "";
  srand(47);
  char alph[] = "ACGT";
  for (int i = 0; i < 50; i++) {
    part1 += alph[rand()%4];
    part2 += alph[rand()%4];
    part3 += alph[rand()%4];
  }

  stringstream ss;
  ss << "@a" << endl;
  ss << part1 << part2 << endl;
  ss << "+" << endl;
  ss << part1 << part2 << endl;
  ss << "@b" << endl;
  ss << part2 << part3 << endl;
  ss << "+" << endl;
  ss << part2 << part3 << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss);

  string genome = "ACGTTT" + part1 + part2 + part3 + "GTCT";
  vector<ReadAlignment> als = rs.GetAlignments(genome);
  ASSERT_EQ(2, als.size());
  EXPECT_EQ(false, als[0].reversed);
  EXPECT_EQ(6, als[0].genome_pos);
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);
  EXPECT_EQ(false, als[1].reversed);
  EXPECT_EQ(56, als[1].genome_pos);
  EXPECT_EQ(1, als[1].read_id);
  EXPECT_EQ(0, als[1].dist);

  genome = "ACGTTT" + ReverseSeq(part1+part2+part3) + "GTCT";
  als = rs.GetAlignments(genome);
  ASSERT_EQ(2, als.size());
  EXPECT_EQ(true, als[0].reversed);
  EXPECT_EQ(56, als[0].genome_pos);
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(0, als[0].dist);
  EXPECT_EQ(true, als[1].reversed);
  EXPECT_EQ(6, als[1].genome_pos);
  EXPECT_EQ(1, als[1].read_id);
  EXPECT_EQ(0, als[1].dist);

  genome = "ACGTTT" + part1 + "ACTGAA" + part2 + "GGAATT" + part3 + "GTCT";
  als = rs.GetAlignments(genome);
  ASSERT_EQ(2, als.size());
  EXPECT_EQ(false, als[0].reversed);
  EXPECT_EQ(6, als[0].genome_pos);
  EXPECT_EQ(0, als[0].read_id);
  EXPECT_EQ(6, als[0].dist);
  EXPECT_EQ(false, als[1].reversed);
  EXPECT_EQ(62, als[1].genome_pos);
  EXPECT_EQ(1, als[1].read_id);
  EXPECT_EQ(6, als[1].dist);
}
