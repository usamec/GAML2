#include "moves.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>

TEST(MovesTest, ForceJoinTest) {
  Node* a = new Node;
  Node* ar = new Node;
  a->id_ = 1;
  a->str_ = "AAAAA";
  ar->str_ = "TTTTT";
  a->rc_ = ar;
  ar->rc_ = a;
  Node* b = new Node;
  Node* br = new Node;
  b->id_ = 2;
  b->str_ = "CCCCC";
  br->str_ = "GGGGG";
  b->rc_ = br;
  br->rc_ = b;
  a->AddNext(b);
  br->AddNext(ar);
  Path p1({a});
  Path p2({b});
  vector<Path> paths({p1, p2});
  vector<Path> out_paths;
  MoveConfig config;
  config.big_node_threshold = 5;
  MakeMove(paths, out_paths, config);
  EXPECT_EQ(1, out_paths.size());
  EXPECT_EQ(a, out_paths[0][0]);
  EXPECT_EQ(b, out_paths[0][1]);
}
