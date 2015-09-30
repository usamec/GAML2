#include "path.h"
#include "graph.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>

TEST(PathTest, BuildPathTest) {
  Node* a = new Node;
  a->id_ = 1;
  Node* b = new Node;
  b->id_ = 2;
  vector<Node*> nodes({a, b});
  vector<Path> paths = BuildPathsFromSingleNodes(nodes);
  ASSERT_EQ(2, paths.size());
  ASSERT_EQ(a, paths[0][0]);
  ASSERT_EQ(b, paths[1][0]);
}

TEST(PathTest, CheckPathTest1) {
  Node* a = new Node;
  a->id_ = 1;
  Node* b = new Node;
  b->id_ = 2;
  a->next_.push_back(b);
  Path p({a, b});
  EXPECT_EQ(true, p.CheckPath());
}

TEST(PathTest, CheckPathTest2) {
  Node* a = new Node;
  a->id_ = 1;
  Node* b = new Node;
  b->id_ = 2;
  Path p({a, b});
  EXPECT_EQ(false, p.CheckPath());
}

TEST(PathTest, TestPathOut) {
  Node* a = new Node;
  a->id_ = 1;
  Node* b = new Node;
  b->id_ = 2;
  Path p({a, b});
  EXPECT_EQ("(1,2)", p.ToDebugString());
}

TEST(PathTest, AppendPathTest) {
  Node* a = new Node;
  a->id_ = 1;
  Node* b = new Node;
  b->id_ = 2;
  Path p({a, b});
  Node* c = new Node;
  c->id_ = 3;
  Node* d = new Node;
  d->id_ = 4;
  Path p2({c, d});
  p.AppendPath(p2);
  ASSERT_EQ(4, p.size());
  EXPECT_EQ(a, p[0]);
  EXPECT_EQ(b, p[1]);
  EXPECT_EQ(c, p[2]);
  EXPECT_EQ(d, p[3]);
}

TEST(PathTest, ReversePathTest1) {
  Node* a = new Node;
  a->id_ = 1;
  Node* b = new Node;
  b->id_ = 2;
  a->rc_ = b;
  b->rc_ = a;
  Node* c = new Node;
  c->id_ = 3;
  Node* d = new Node;
  d->id_ = 4;
  c->rc_ = d;
  d->rc_ = c;
  Path p({a, c});
  p.Reverse();
  EXPECT_EQ(4, p[0]->id_);
  EXPECT_EQ(2, p[1]->id_);
  Path p2 = p.GetReverse();
  EXPECT_EQ(1, p2[0]->id_);
  EXPECT_EQ(3, p2[1]->id_);
}

TEST(PathTest, AppendPathWithGapTest) {
  Node* a = new Node;
  Node* ar = new Node;
  a->id_ = 1;
  a->rc_ = ar;
  ar->rc_ = a;
  Node* b = new Node;
  Node* br = new Node;
  b->id_ = 2;
  b->rc_ = br;
  br->rc_ = b;
  a->AddNext(b);
  br->AddNext(ar);
  Path p({a, b});
  Node* c = new Node;
  Node* cr = new Node;
  c->id_ = 3;
  c->rc_ = cr;
  cr->rc_ = c;
  Node* d = new Node;
  Node* dr = new Node;
  d->id_ = 4;
  d->rc_ = dr;
  dr->rc_ = d;
  c->AddNext(d);
  dr->AddNext(cr);
  Path p2({c, d});
  p.AppendPathWithGap(p2, 47);
  ASSERT_EQ(5, p.size());
  EXPECT_EQ(a, p[0]);
  EXPECT_EQ(b, p[1]);
  EXPECT_EQ(c, p[3]);
  EXPECT_EQ(d, p[4]);
  EXPECT_EQ(true, p.CheckPath());
  EXPECT_EQ(47, p[2]->GapLength());
  EXPECT_EQ(true, p[2]->IsGap());
  p.Reverse();
  ASSERT_EQ(5, p.size());
  EXPECT_EQ(dr, p[0]);
  EXPECT_EQ(cr, p[1]);
  EXPECT_EQ(br, p[3]);
  EXPECT_EQ(ar, p[4]);
  EXPECT_EQ(true, p.CheckPath());
  EXPECT_EQ(47, p[2]->GapLength());
  EXPECT_EQ(true, p[2]->IsGap());
}

TEST(PathTest, PathToStringTest) {
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
  string str_out1 = p1.ToString(false);
  EXPECT_EQ(121, str_out1.size());
  EXPECT_EQ(g->nodes_[0]->str_, str_out1);
  string str_out2 = p1.ToString(true);
  EXPECT_EQ(161, str_out2.size());
  EXPECT_EQ("TTTCATAGAAAGCATTTTGTTGTTCTTTGTTGAATTTGTT"
            "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA",
            str_out2);

  Path p2({g->nodes_[0], g->nodes_[2]});
  string str_out3 = p2.ToString(false);
  EXPECT_EQ(g->nodes_[0]->str_ + g->nodes_[2]->str_, str_out3);
  string str_out4 = p2.ToString(true);
  EXPECT_EQ("TTTCATAGAAAGCATTTTGTTGTTCTTTGTTGAATTTGTT"
            "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTAAGAC",
            str_out4);
}

TEST(PathTest, PathToStringTest2) {
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
  Path p2({g->nodes_[2]});
  p1.AppendPathWithGap(p2, 50);
  string str_out1 = p1.ToString(false);
  EXPECT_EQ(121+50+4, str_out1.size());
  EXPECT_EQ(g->nodes_[0]->str_ + string(50, 'N') + g->nodes_[2]->str_, str_out1);
  string str_out2 = p1.ToString(true);
  EXPECT_EQ(161+50+4, str_out2.size());
  EXPECT_EQ("TTTCATAGAAAGCATTTTGTTGTTCTTTGTTGAATTTGTT"
            "GTCAGCTTTTGGTGCTTGAGCATCATTTAGCTTTTTAGCTTCTGCTAAAAGGTTAGCGCTTTGGCTTGGGTCATCTTTTAGGCTTTGGATGAAACCATTGCGTTGTTCTTCGTTTAAGTTA"
            + string(50, 'N') + "AGAC",
            str_out2);
}
