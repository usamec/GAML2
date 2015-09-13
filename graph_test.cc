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
