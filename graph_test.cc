#include "graph.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>

TEST(GraphTest, VelvetNumToMyNumTest) {
  EXPECT_EQ(0, VelvetNumToMyNum(1));
  EXPECT_EQ(1, VelvetNumToMyNum(-1));
  EXPECT_EQ(2, VelvetNumToMyNum(2));
  EXPECT_EQ(3, VelvetNumToMyNum(-2));
}

TEST(GraphTest, LoadArcTest) {
  stringstream ss;
  ss << "ARC\t3\t-2\t47\n";
  int a, b;
  tie(a, b) = LoadArc(ss);
  EXPECT_EQ(4, a);
  EXPECT_EQ(3, b);
}

TEST(GraphTest, LoadArcTest2) {
  stringstream ss;
  ss << "SS\t3\t-2\t47\n";
  int a, b;
  tie(a, b) = LoadArc(ss);
  EXPECT_EQ(-1, a);
  EXPECT_EQ(-1, b);
}

TEST(GraphTest, LoadArcTest3) {
  stringstream ss;
  int a, b;
  tie(a, b) = LoadArc(ss);
  EXPECT_EQ(-1, a);
  EXPECT_EQ(-1, b);
}

TEST(GraphTest, LoadGraphTest) {
  stringstream ss;
  ss << "2\t1000\t41\t1\n";
  ss << "NODE\t1\t4\t0\t0\n";
  ss << "AAAC\n";
  ss << "TTCC\n";
  ss << "NODE\t2\t4\t0\t0\n";
  ss << "AGAC\n";
  ss << "TGCC\n";
  ss << "ARC\t1\t-2\t44\n";
  Graph *g = LoadGraph(ss);

  EXPECT_EQ(4, g->nodes_.size());
  EXPECT_EQ(41, g->k_);
  EXPECT_EQ(g->nodes_[1], g->nodes_[0]->rc_);
  EXPECT_EQ(g->nodes_[0], g->nodes_[1]->rc_);
  EXPECT_EQ(g->nodes_[3], g->nodes_[2]->rc_);
  EXPECT_EQ(g->nodes_[2], g->nodes_[3]->rc_);
  EXPECT_EQ("AAAC", g->nodes_[0]->str_);
  EXPECT_EQ("TTCC", g->nodes_[1]->str_);
  EXPECT_EQ("AGAC", g->nodes_[2]->str_);
  EXPECT_EQ("TGCC", g->nodes_[3]->str_);
  ASSERT_EQ(1, g->nodes_[0]->next_.size());
  ASSERT_EQ(0, g->nodes_[1]->next_.size());
  ASSERT_EQ(1, g->nodes_[2]->next_.size());
  ASSERT_EQ(0, g->nodes_[3]->next_.size());
  ASSERT_EQ(g->nodes_[3], g->nodes_[0]->next_[0]);
  ASSERT_EQ(g->nodes_[1], g->nodes_[2]->next_[0]);

  EXPECT_EQ(g, g->nodes_[0]->graph_);
  EXPECT_EQ(g, g->nodes_[1]->graph_);
  EXPECT_EQ(g, g->nodes_[2]->graph_);
  EXPECT_EQ(g, g->nodes_[3]->graph_);
}
