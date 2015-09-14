#include "util.h"
#include <gtest/gtest.h>

TEST(ReverseTest, ReverseTest) {
  EXPECT_EQ(string("TTACG"), ReverseSeq("CGTAA"));
}
