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
