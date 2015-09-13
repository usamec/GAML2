#include "node.h"
#include <gtest/gtest.h>
#include <sstream>
#include <tuple>

TEST(NodeTest, LoadNodeTest) {
  stringstream ss;
  ss << "NODE\t1\t4\t99998\t99998\t0\t0" << endl;
  ss << "AACG" << endl;
  ss << "TTCC" << endl;

  Node *f, *b;
  tie(f, b) = LoadNode(ss, 5);
  EXPECT_EQ(10, f->id_);
  EXPECT_EQ(11, b->id_);
  EXPECT_EQ("AACG", f->str_);
  EXPECT_EQ("TTCC", b->str_);
  EXPECT_EQ(b, f->rc_);
  EXPECT_EQ(f, b->rc_);
}
