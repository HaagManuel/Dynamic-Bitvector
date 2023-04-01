#include <gtest/gtest.h>
#include <vector>

#include "src/static_small_bv.h"
#include "src/naive_bv.h"
#include "src/random_generator.h"


//test
#include "src/dynamic_small_bv.h"
//test

// #define STATIC_BV StaticSmallBV
#define STATIC_BV DynamicSmallBV


TEST(TestStaticBV, SizeTest) {
  STATIC_BV<4> bv;
  EXPECT_EQ(bv.size(), 0);
  EXPECT_EQ(bv.max_size(), 4 * 64);
  bv.insert(0, true);
  EXPECT_EQ(bv.size(), 1);
  bv.insert(1, true);
  EXPECT_EQ(bv.size(), 2);
  bv.insert(2, true);
  EXPECT_EQ(bv.size(), 3);
}

TEST(TestStaticBV, InsertBlock) {
  STATIC_BV<4> bv;
  uint32_t N = 130;
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(i, true);
  }
  EXPECT_EQ(bv.size(), N);
  for(uint32_t i = 0; i < N; i++) {
    EXPECT_EQ(bv.access(i), true);
  }
  // for(uint32_t i = N; i < bv.max_size(); i++) {
  //   EXPECT_EQ(bv.access(i), false);
  // }
}

TEST(TestStaticBV, InsertMod2) {
  STATIC_BV<4> bv;
  uint32_t N = bv.max_size() - 4;
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(i, i % 2 == 0);
  }
  EXPECT_EQ(bv.size(), N);
  for(uint32_t i = 0; i < N; i++) {
    EXPECT_EQ(bv.access(i), i % 2 == 0);
  }
  // for(uint32_t i = N; i < bv.max_size(); i++) {
  //   EXPECT_EQ(bv.access(i), false);
  // }
}

TEST(TestStaticBV, InsertMod3Front) {
  STATIC_BV<4> bv;
  uint32_t N = 100; // % 3 == 1 -> alignment
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(0, i % 3 == 0);
  }
  ASSERT_EQ(bv.size(), N);
  for(uint32_t i = 0; i < N; i++) {
    ASSERT_EQ(bv.access(i), i % 3 == 0);
  }
  // for(uint32_t i = N; i < bv.max_size(); i++) {
  //   ASSERT_EQ(bv.access(i), false);
  // }
}

TEST(TestStaticBV, InsertTest3) {
  STATIC_BV<3> bv;
  uint32_t N = 130;
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(i, i % 5 == 0);
  }
  EXPECT_EQ(bv.size(), N);
  for(uint32_t i = 0; i < N; i++) {
    EXPECT_EQ(bv.access(i), i % 5 == 0);
  }
  // for(uint32_t i = N; i < bv.max_size(); i++) {
  //   EXPECT_EQ(bv.access(i), false);
  // }
}

TEST(TestStaticBV, InsertMid) {
  STATIC_BV<3> bv;
  uint32_t N = 10;
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(i, true);
  }
  bv.insert(5, false);
  for(uint32_t i = 0; i < N + 1; i++) {
    EXPECT_EQ(bv.access(i), i != 5);
  }
}

TEST(TestStaticBV, RankSelectTest1) {
  STATIC_BV<2> bv;
  uint32_t N = 70;
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(i, true);
  }
  for(uint32_t i = 0; i < N; i++) {
    ASSERT_EQ(bv.rank(i, true), i);
    ASSERT_EQ(bv.select(i + 1, true), i);
    ASSERT_EQ(bv.rank(i, false), 0);
    ASSERT_EQ(bv.select(i + 1, false), bv.max_size());
  }
  ASSERT_EQ(bv.select(N + 1, true), bv.max_size());
  ASSERT_EQ(bv.select(N + 1, false), bv.max_size());
}

TEST(TestStaticBV, RandomOp) {
  STATIC_BV<32> bv;
  NaiveBV bv2;
  //random is always same sequence
  uint32_t N = bv.max_size();
  for(uint32_t i = 0; i < N; i++) {
    uint32_t j = random_gen::random_index(i + 1);
    bool b = random_gen::random_bool();
    bv.insert(j, b);
    bv2.insert(j, b);
    ASSERT_EQ(bv.size(), bv2.size());
  }
  for(uint32_t i = 0; i < N; i++) {
    uint32_t j = random_gen::random_index(N);
    uint32_t k = random_gen::random_index(N);
    uint32_t l = random_gen::random_index(N);
    bool b = random_gen::random_bool();
    ASSERT_EQ(bv.access(j), bv2.access(j));
    ASSERT_EQ(bv.select(j + 1, b), bv2.select(j + 1, b));
    ASSERT_EQ(bv.rank(j , b), bv2.rank(j, b));
    bv.flip(j);
    bv2.flip(j);
    ASSERT_EQ(bv.access(j), bv2.access(j));
    bv.remove(k);
    bv2.remove(k);
    bv.insert(l, b);
    bv2.insert(l, b);
    ASSERT_EQ(bv.access(l), bv2.access(l));

  }
}

TEST(TestStaticBV, RemoveFrontMod2) {
  STATIC_BV<3> bv;
  uint32_t N = 80;
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(i, i % 2 == 0);
  }
  ASSERT_EQ(bv.size(), N);
  for(uint32_t i = 0; i < N; i++) {
    bv.remove(0);
    ASSERT_EQ(bv.size(), N - (i + 1));
    for(uint32_t j = 0; j < N - (i + 1); j++) {
      ASSERT_EQ(bv.access(j), (i + j) % 2 == 1);
    }
  }
}

TEST(TestStaticBV, SplitTest) {
  uint32_t N = 4 * 64;
  ASSERT_EQ(N % 2, 0);
  std::vector<bool> seq = random_gen::random_bit_string(N);
  STATIC_BV<4> bv;
  STATIC_BV<4> bv2;
  for(uint32_t i = 0; i < N; i++) {
    bv.insert(i, seq[i]);
  }
  ASSERT_EQ(bv.size(), N);
  ASSERT_EQ(bv2.size(), 0);
  bv.split(bv2);
  ASSERT_EQ(bv.size(), N / 2);
  ASSERT_EQ(bv2.size(), N / 2);
  for(uint32_t i = 0; i < N / 2; i++) {
    ASSERT_EQ(bv.access(i), seq[i]);
    ASSERT_EQ(bv2.access(i), seq[i +  N / 2]);
  }
}

TEST(TestStaticBV, DISABLED_AppendTest) {
  uint32_t N = 64;
  uint32_t off = 4;
  std::vector<bool> seq = random_gen::random_bit_string(8 * N);
  STATIC_BV<8> bv;
  STATIC_BV<8> bv2;
  for(uint32_t i = 0; i < 3 * N; i++) {
    bv.insert(i, seq[i]);
  }
  for(uint32_t i = 3 * N; i < 7 * N + off; i++) {
    bv2.insert(i - 3 * N, seq[i]);
  }
  ASSERT_EQ(bv.size(), 3 * N);
  ASSERT_EQ(bv2.size(), 7 * N + off - 3 * N);
  ASSERT_EQ(bv.is_aligned(), true);
  ASSERT_EQ(bv2.is_aligned(), false);
  bv.append_bv(bv2);
  ASSERT_EQ(bv.size(), 7 * N + off);
  ASSERT_EQ(bv2.size(), 7 * N + off - 3 * N);
 for(uint32_t i = 0; i < 7 * N + off; i++) {
    ASSERT_EQ(bv.access(i), seq[i]);
  }
}

//DISABLED_ -> skip
TEST(TestStaticBV, DISABLED_BalanceLeftRight) {
  uint32_t N = 64;
  uint32_t off = 4;
  uint32_t left = 6 * N, right = N + off;
  std::vector<bool> seq = random_gen::random_bit_string(8 * N);
  STATIC_BV<8> bv;
  STATIC_BV<8> bv2;
  for(uint32_t i = 0; i < left; i++) {
    bv.insert(i, seq[i]);
  }
  for(uint32_t i = left; i < left + right; i++) {
    bv2.insert(i - left, seq[i]);
  }
  ASSERT_EQ(bv.size(), left);
  ASSERT_EQ(bv2.size(), right);
  ASSERT_EQ(bv.is_aligned(), true);
  ASSERT_EQ(bv2.is_aligned(), false);
  ASSERT_EQ(bv.last_used_block(), 5);
  ASSERT_EQ(bv2.last_used_block(), 1);
  
  STATIC_BV<8>::balance_left_to_right(bv, bv2);

  ASSERT_EQ(bv.size(), left - 2 * N);
  ASSERT_EQ(bv2.size(), right + 2 * N);
 for(uint32_t i = 0; i < left - 2 * N; i++) {
    ASSERT_EQ(bv.access(i), seq[i]);
  }
  for(uint32_t i = left - 2 * N; i < left + right; i++) {
    ASSERT_EQ(bv2.access(i - (left - 2 * N)), seq[i]);
  }
}

TEST(TestStaticBV, DISABLED_BalanceRightLeft1) {
  uint32_t N = 64;
  uint32_t off = 4;
  uint32_t left = N, right = 6 * N + off;
  std::vector<bool> seq = random_gen::random_bit_string(8 * N);
  STATIC_BV<8> bv;
  STATIC_BV<8> bv2;
  for(uint32_t i = 0; i < left; i++) {
    bv.insert(i, seq[i]);
  }
  for(uint32_t i = left; i < left + right; i++) {
    bv2.insert(i - left, seq[i]);
  }
  ASSERT_EQ(left < right, true);
  ASSERT_EQ(bv.size(), left);
  ASSERT_EQ(bv2.size(), right);
  ASSERT_EQ(bv.is_aligned(), true);
  ASSERT_EQ(bv2.is_aligned(), false);
  ASSERT_EQ(bv.last_used_block(), 0);
  ASSERT_EQ(bv2.last_used_block(), 6);
  
  STATIC_BV<8>::balance_right_to_left(bv, bv2);

  ASSERT_EQ(bv.size(), left + 3 * N);
  ASSERT_EQ(bv2.size(), right - 3 * N);
 for(uint32_t i = 0; i < left + 3 * N; i++) {
    ASSERT_EQ(bv.access(i), seq[i]);
  }
  for(uint32_t i = left + 3 * N; i < left + right; i++) {
    ASSERT_EQ(bv2.access(i - (left + 3 * N)), seq[i]);
  }
}

TEST(TestStaticBV, DISABLED_BalanceRightLeft2) {
  uint32_t N = 64;
  uint32_t off = 4;
  uint32_t left = N + off, right = 6 * N + off;
  std::vector<bool> seq = random_gen::random_bit_string(8 * N);
  STATIC_BV<8> bv;
  STATIC_BV<8> bv2;
  for(uint32_t i = 0; i < left; i++) {
    bv.insert(i, seq[i]);
  }
  for(uint32_t i = left; i < left + right; i++) {
    bv2.insert(i - left, seq[i]);
  }
  ASSERT_EQ(left < right, true);
  ASSERT_EQ(bv.size(), left);
  ASSERT_EQ(bv2.size(), right);
  ASSERT_EQ(bv.is_aligned(), false);
  ASSERT_EQ(bv2.is_aligned(), false);
  ASSERT_EQ(bv.last_used_block(), 1);
  ASSERT_EQ(bv2.last_used_block(), 6);
  
  STATIC_BV<8>::balance_right_to_left(bv, bv2);

  ASSERT_EQ(bv.size(), left + 3 * N - off);
  ASSERT_EQ(bv2.size(), right - 3 * N + off);
 for(uint32_t i = 0; i < left + 3 * N - off; i++) {
    ASSERT_EQ(bv.access(i), seq[i]);
  }
  for(uint32_t i = left + 3 * N - off; i < left + right; i++) {
    ASSERT_EQ(bv2.access(i - (left + 3 * N - off)), seq[i]);
  }
}