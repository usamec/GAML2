#include "read_set.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>
#include <set>
#include <algorithm>
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

  SingleShortReadSet<> rs;
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

  SingleShortReadSet<> rs;
  rs.LoadReadSet(ss);

  CandidateReadPosition candidate(0, 0, 0);
  SingleReadAlignment al;
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

  SingleShortReadSet<> rs;
  rs.LoadReadSet(ss);

  vector<SingleReadAlignment> als = rs.GetAlignments("AAAACCCCTTTTGGGG");
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

  SingleShortReadSet<> rs;
  rs.LoadReadSet(ss);

  string genome = "ACGTTT" + part1 + part2 + "GTCT";
  vector<SingleReadAlignment> als = rs.GetAlignments(genome);
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

  SingleShortReadSet<> rs;
  rs.LoadReadSet(ss);

  string genome = "ACGTTT" + part1 + part2 + part3 + "GTCT";
  vector<SingleReadAlignment> als = rs.GetAlignments(genome);
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

char GetRandomBase(char otherThan = 'X') {
  char c = 'A';
  do {
    switch (rand() % 4) {
      case 0: c = 'A'; break;
      case 1: c = 'C'; break;
      case 2: c = 'G'; break;
      case 3: c = 'T'; break;
    }
  } while (c == otherThan);
  return c;
}

string RemapReverse(const string& sequence) {
  stringstream resultStream;
  for (char c : sequence) {
    switch (c) {
      case 'A': resultStream << 'T'; break;
      case 'C': resultStream << 'G'; break;
      case 'G': resultStream << 'C'; break;
      case 'T': resultStream << 'A'; break;
    }
  }
  string result = resultStream.str();
  reverse(result.begin(), result.end());
  return result;
}

TEST(ReadSetTest, PacBioTest) {
  srand(47);
  
  const int genomeLength = 10000;
  const int readLength = 1000;
  const int numReads = 10;
  const int readPadding = 400;
  const int maxAlignmentEndDifference = 100;
  const int maxDistanceDifference = 40;
  
  // error rates (n out of 100)
  const int insertRate = 4;
  const int deleteRate = 2;
  const int substituteRate = 1;
  
  // generate genome
  stringstream genomeStream;
  for (int i = 0; i < genomeLength; i++) {
    genomeStream << GetRandomBase();
  }
  string genome = genomeStream.str();
  
  vector<string> reads; // reads strings
  vector<int> startPositions; // starting position of valid part of read in the genome
  vector<bool> areReversed;
  vector<int> numErrors;
  
  stringstream fastqStream;
  
  for (int j = 0; j < numReads; j++) {
    int readPos = rand() % (genomeLength - readLength + 1);
    
    stringstream readStream;
    // add padding
    for (int i = 0; i < readPadding; i++) {
      readStream << GetRandomBase();
    }
    int errors = 0;
    // create noisy read
    for (int i = readPos; i < readPos + readLength; i++) {
      if (rand() % 100 < insertRate) {
        readStream << GetRandomBase();
        i--;
        errors++;
      } else if (rand() % 100 < deleteRate) {
        // nothing to do here
        errors++;
      } else if (rand() % 100 < substituteRate) {
        readStream << GetRandomBase(genome[i]);
        errors++;
      } else {
        readStream << genome[i];
      }
    }
    // add padding
    for (int i = 0; i < readPadding; i++) {
      readStream << GetRandomBase();
    }
    
    string read = readStream.str();
    
    // maybe reverse
    bool reversed = (rand() % 2) == 0;
    if (reversed) {
      read = RemapReverse(read);
    }
    
    reads.push_back(read);
    startPositions.push_back(readPos);
    areReversed.push_back(reversed);
    numErrors.push_back(errors);
    
    fastqStream << "@read" << j << endl;
    fastqStream << read << endl;
    fastqStream << "+" << endl;
    fastqStream << read << endl;
  }
  
  ReadSetPacBio<StandardReadIndex> rs;
  rs.SetParameters(0.90, {0.25, 0.25, 0.25, 0.25}, 500);
  rs.LoadReadSet(fastqStream);
  
  vector<ReadAlignmentPacBio> alignments = rs.GetAlignments(genome);
  
  set<int> alignedReadsIds;
  for (unsigned int i = 0; i < alignments.size(); i++) {
    ReadAlignmentPacBio &al = alignments[i];
    alignedReadsIds.insert(al.read_id);
    
    // check position withing genome
    EXPECT_NEAR(al.genome_first, startPositions[al.read_id], maxAlignmentEndDifference);
    EXPECT_NEAR(al.genome_last, startPositions[al.read_id] + readLength, maxAlignmentEndDifference);
    
    // check position within read
    EXPECT_NEAR(al.read_first, readPadding, maxAlignmentEndDifference);
    EXPECT_NEAR(al.read_last, (int)reads[al.read_id].length() - readPadding, maxAlignmentEndDifference);
    
    // check orientation
    EXPECT_EQ(al.reversed, areReversed[al.read_id]);
    
    // check difference
    EXPECT_NEAR(al.dist, numErrors[al.read_id], maxDistanceDifference);
  }
  
  // all reads must be found
  EXPECT_EQ(numReads, alignedReadsIds.size());
}
