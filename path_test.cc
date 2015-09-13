#include "path.h"
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
