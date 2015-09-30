#include "read_probability_calculator.h"
#include "graph.h"
#include "util.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>


TEST(SingleReadProbabilityCalculatorTest, Test1) {
  stringstream ss;
  ss << "2\t1000\t41\t1\n";
  ss << "NODE\t1\t121\t0\t0\n";
  ss << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA\n";
  ss << "CTAAAAGATGACCCAAGCCAAAGCGCTAACCTTTTAGCAGAAGCTAAAAAGCTAAATGATGCTCAAGCACCAAAAGCTGACAACAAATTCAACAAAGAACAACAAAATGCTTTCTATGAAA\n";
  ss << "NODE\t2\t4\t0\t0\n";
  ss << "AGAC\n";
  ss << "TGCC\n";
  ss << "ARC\t1\t-2\t44\n";
  Graph *g = LoadGraph(ss);

  Path p1({g->nodes_[0]});

  stringstream ss2;
  ss2 << "@a" << endl;
  ss2 << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA" << endl;
  ss2 << "+" << endl;
  ss2 << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA" << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss2);

  SingleReadProbabilityCalculator rp(&rs, 0.01, -10, -0.7, 0, 0);
  ProbabilityChange prob_change;
  double pr1 = rp.GetPathsProbability(vector<Path>({p1}), prob_change);
  EXPECT_FLOAT_EQ(-6.29749500325813737913, pr1);
}

TEST(SingleReadProbabilityCalculatorTest, Test2) {
  stringstream ss;
  ss << "2\t1000\t41\t1\n";
  ss << "NODE\t1\t121\t0\t0\n";
  ss << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA\n";
  ss << "CTAAAAGATGACCCAAGCCAAAGCGCTAACCTTTTAGCAGAAGCTAAAAAGCTAAATGATGCTCAAGCACCAAAAGCTGACAACAAATTCAACAAAGAACAACAAAATGCTTTCTATGAAA\n";
  ss << "NODE\t2\t4\t0\t0\n";
  ss << "AGAC\n";
  ss << "TGCC\n";
  ss << "ARC\t1\t-2\t44\n";
  Graph *g = LoadGraph(ss);

  Path p1({g->nodes_[0]});

  stringstream ss2;
  ss2 << "@a" << endl;
  ss2 << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
  ss2 << "+" << endl;
  ss2 << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA" << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss2);

  SingleReadProbabilityCalculator rp(&rs, 0.01, -10, -0.7, 0, 0);
  ProbabilityChange prob_change;
  double pr1 = rp.GetPathsProbability(vector<Path>({p1}), prob_change);
  EXPECT_FLOAT_EQ(-99.781403, pr1);
}

TEST(SingleReadProbabilityCalculatorTest, Test3) {
  string part1 = "";
  char alph[] = "ACGT";
  for (int i = 0; i < 50; i++) {
    part1 += alph[rand()%4];
  }
  stringstream ss;
  ss << "2\t1000\t41\t1\n";
  ss << "NODE\t1\t180\t0\t0\n";
  ss << part1 + string(40, 'C') + part1 + string(40, 'A') << endl;
  ss << ReverseSeq(part1) + string(40, 'G') + ReverseSeq(part1) + string(40, 'A') << endl;
  ss << "NODE\t2\t4\t0\t0\n";
  ss << "AGAC\n";
  ss << "TGCC\n";
  ss << "ARC\t1\t-2\t44\n";
  Graph *g = LoadGraph(ss);

  Path p1({g->nodes_[0]});

  stringstream ss2;
  ss2 << "@a" << endl;
  ss2 << part1 << endl;
  ss2 << "+" << endl;
  ss2 << part1 << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss2);

  SingleReadProbabilityCalculator rp(&rs, 0.01, -10, -0.7, 0, 0);
  ProbabilityChange prob_change;
  double pr1 = rp.GetPathsProbability(vector<Path>({p1}), prob_change);
  EXPECT_FLOAT_EQ(-5.202997, pr1);
}

TEST(SingleReadProbabilityCalculatorTest, Test4) {
  stringstream ss;
  ss << "2\t1000\t41\t1\n";
  ss << "NODE\t1\t121\t0\t0\n";
  ss << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA\n";
  ss << "CTAAAAGATGACCCAAGCCAAAGCGCTAACCTTTTAGCAGAAGCTAAAAAGCTAAATGATGCTCAAGCACCAAAAGCTGACAACAAATTCAACAAAGAACAACAAAATGCTTTCTATGAAA\n";
  ss << "NODE\t2\t4\t0\t0\n";
  ss << "AGAC\n";
  ss << "TGCC\n";
  ss << "ARC\t1\t-2\t44\n";
  Graph *g = LoadGraph(ss);

  Path p1({g->nodes_[0]});

  stringstream ss2;
  ss2 << "@a" << endl;
  ss2 << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA" << endl;
  ss2 << "+" << endl;
  ss2 << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA" << endl;
  ss2 << "@b" << endl;
  ss2 << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA" << endl;
  ss2 << "+" << endl;
  ss2 << "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA" << endl;

  ReadSet<> rs;
  rs.LoadReadSet(ss2);

  SingleReadProbabilityCalculator rp(&rs, 0.01, -10, -0.7, 0, 0);
  ProbabilityChange prob_change;
  double pr1 = rp.GetPathsProbability(vector<Path>({p1}), prob_change);
  EXPECT_FLOAT_EQ(-6.29749500325813737913, pr1);
}
