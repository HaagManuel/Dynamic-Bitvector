#include <gtest/gtest.h>

#include "src/naive_bv.h"

TEST(TestNaiveBV, SizeCheck) {
  NaiveBV bv;
  bv.insert(0, true);
  EXPECT_EQ(bv.size(), 1);
  bv.insert(0, true);
  EXPECT_EQ(bv.size(), 2);
  bv.remove(0);
  EXPECT_EQ(bv.size(), 1);
  bv.remove(0);
  EXPECT_EQ(bv.size(), 0);
}

TEST(TestNaiveBV, InsertCheck) {
  NaiveBV bv;
  bv.insert(0, true);
  bv.insert(1, true);
  bv.insert(2, false);
  bv.insert(3, false);
  EXPECT_EQ(bv.size(), 4);
  EXPECT_EQ(bv.access(0), true);
  EXPECT_EQ(bv.access(1), true);
  EXPECT_EQ(bv.access(2), false);
  EXPECT_EQ(bv.access(3), false);
}

TEST(TestNaiveBV, RemoveCheck) {
  NaiveBV bv;
  bv.insert(0, true);
  bv.insert(1, true);
  bv.insert(2, false);
  bv.insert(3, false);
  bv.remove(2);
  EXPECT_EQ(bv.size(), 3);
  EXPECT_EQ(bv.access(0), true);
  EXPECT_EQ(bv.access(1), true);
  EXPECT_EQ(bv.access(2), false);
  bv.remove(0);
  EXPECT_EQ(bv.size(), 2);
  EXPECT_EQ(bv.access(0), true);
  EXPECT_EQ(bv.access(1), false);
}

TEST(TestNaiveBV, RankSelectCheck) {
  NaiveBV bv;
  bv.insert(0, true);
  bv.insert(1, false);
  bv.insert(2, true);
  bv.insert(3, false);
  bv.insert(4, true);
  EXPECT_EQ(bv.size(), 5);

  EXPECT_EQ(bv.rank(0, true), 0);
  EXPECT_EQ(bv.rank(1, true), 1);
  EXPECT_EQ(bv.rank(2, true), 1);
  EXPECT_EQ(bv.rank(3, true), 2);
  EXPECT_EQ(bv.rank(4, true), 2);

  EXPECT_EQ(bv.rank(0, false), 0);
  EXPECT_EQ(bv.rank(1, false), 0);
  EXPECT_EQ(bv.rank(2, false), 1);
  EXPECT_EQ(bv.rank(3, false), 1);
  EXPECT_EQ(bv.rank(4, false), 2);

  EXPECT_EQ(bv.select(1, true), 0);
  EXPECT_EQ(bv.select(2, true), 2);
  EXPECT_EQ(bv.select(3, true), 4);
  EXPECT_EQ(bv.select(4, true), bv.size());

  EXPECT_EQ(bv.select(1, false), 1);
  EXPECT_EQ(bv.select(2, false), 3);
  EXPECT_EQ(bv.select(3, false), bv.size());
}