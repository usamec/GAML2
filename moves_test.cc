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
  MakeMove(paths, out_paths, config, false);
  EXPECT_EQ(1, out_paths.size());
  EXPECT_TRUE(
      (a == out_paths[0][0] && b == out_paths[0][1]) ||
      (br == out_paths[0][0] && ar == out_paths[0][1]));
}

TEST(MovesTest, ForceBreakTest) {
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
  b->str_ = "C";
  br->str_ = "G";
  b->rc_ = br;
  br->rc_ = b;
  a->AddNext(b);
  br->AddNext(ar);
  Node* c = new Node;
  Node* cr = new Node;
  c->id_ = 3;
  c->str_ = "TCCCC";
  cr->str_ = "GGGGA";
  c->rc_ = cr;
  cr->rc_ = c;
  b->AddNext(c);
  cr->AddNext(br);
  Path p({a, b, c});
  vector<Path> paths({p});
  vector<Path> out_paths;
  MoveConfig config;
  config.big_node_threshold = 5;
  MakeMove(paths, out_paths, config, false);
  EXPECT_EQ(2, out_paths.size());
  EXPECT_EQ(1, out_paths[0].size());
  EXPECT_EQ(1, out_paths[1].size());
  EXPECT_EQ(a, out_paths[0][0]);
  EXPECT_EQ(c, out_paths[1][0]);
}
