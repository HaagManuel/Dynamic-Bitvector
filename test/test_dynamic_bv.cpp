#include <gtest/gtest.h>
#include <vector>

#include "src/dynamic_bv.h"
#include "src/naive_bv.h"
#include "src/random_generator.h"

#include "src/oop/dynamic_bv_oop.h"
#include "src/variant/dynamic_bv_variant.h"

#define BV DynamicBV
// #define BV DynamicBVOop
// #define BV DynamicBVVariant


TEST(TestDynamicBV, OnlyRoot) {
  BV<2> bv;
  uint32_t N = 127;
  for(uint32_t i = 0; i < N; i++) {
      bv.insert(i, i % 3 == 0);
      // ASSERT_EQ(bv.size(), i + 1);
  }
  for(uint32_t i = 0; i < N; i++) {
      ASSERT_EQ(bv.access(i), i % 3 == 0);
  }
}

TEST(TestDynamicBV, OneSplit) {
  BV<2> bv;
  uint32_t N = 200;
  for(uint32_t i = 0; i < N; i++) {
      bv.insert(i, i % 5 == 0);
      // ASSERT_EQ(bv.size(), i + 1);
  }
  for(uint32_t i = 0; i < N; i++) {
      ASSERT_EQ(bv.access(i), i % 5 == 0);
  }
}

TEST(TestDynamicBV, TestBalancingRightRight) {
  BV<2> bv;
  uint32_t N = 400;
  for(uint32_t i = 0; i < N; i++) {
      bv.insert(i, true);
      ASSERT_EQ(bv.check_invariant(), true);
  }
}

TEST(TestDynamicBV, TestBalancingLeftLeft) {
  BV<2> bv;
  uint32_t N = 400;
  for(uint32_t i = 0; i < N; i++) {
      bv.insert(0, true);
      ASSERT_EQ(bv.check_invariant(), true);
  }
}

TEST(TestDynamicBV, TestBalancingMid) {
  BV<2> bv;
  uint32_t N = 1000;
  for(uint32_t i = 0; i < N; i++) {
    // std::cout << V(i) << "\n";
    bv.insert(i / 2, true);
    ASSERT_EQ(bv.check_invariant(), true);
  }
}

TEST(TestDynamicBV, RandomFlips) {
  BV<2> bv;
  NaiveBV bv2;
  // random is always same sequence
  uint32_t N = 4000;
  uint32_t count[2] = {0,0};
  for(uint32_t i = 0; i < N; i++) {
    uint32_t j = random_gen::random_index(i + 1);
    bool b = random_gen::random_bool();
    bv.insert(j, b);
    bv2.insert(j, b);
    // ASSERT_EQ(bv.size(), bv2.size());
    count[b]++;
  }
  for(uint32_t i = 0; i < N; i++) {
    ASSERT_EQ(bv.access(i), bv2.access(i));
  }
  for(uint32_t i = 0; i < N; i++) {
    uint32_t j = random_gen::random_index(N);
    count[0] += bv2.access(j) == 0? -1 : 1;
    count[1] += bv2.access(j) == 0? 1 : -1;
    uint32_t k = random_gen::random_index(N);
    uint32_t l1 = random_gen::random_index(count[0]) + 1;
    uint32_t l2 = random_gen::random_index(count[1]) + 1;
    bv.flip(j);
    bv2.flip(j);
    ASSERT_EQ(bv.access(j), bv2.access(j));
    ASSERT_EQ(bv.rank(k, true), bv2.rank(k, true));
    ASSERT_EQ(bv.rank(k, false), bv2.rank(k, false));
    ASSERT_EQ(bv.select(l1, false), bv2.select(l1, false));
    ASSERT_EQ(bv.select(l2, true), bv2.select(l2, true));
  }
  // ASSERT_EQ(bv.select(count[0] + 2, false), bv.size());
  // ASSERT_EQ(bv.select(count[1] + 2, true), bv.size());
}

TEST(TestDynamicBV, RandomWithoutDelete) {
  BV<2> bv;
  NaiveBV bv2;
  // random is always same sequence
  uint32_t N = 4000;
  uint32_t count[2] = {0,0};
  for(uint32_t i = 0; i < N; i++) {
    uint32_t j = random_gen::random_index(i + 1);
    bool b = random_gen::random_bool();
    bv.insert(j, b);
    bv2.insert(j, b);
    // ASSERT_EQ(bv.size(), bv2.size());
    ASSERT_EQ(bv.check_invariant(), true);
    count[b]++;
  }
  // bv.pr_tree();
  for(uint32_t i = 0; i < N; i++) {
    ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    ASSERT_EQ(bv.rank(i, true), bv2.rank(i, true)) << V(i);
    ASSERT_EQ(bv.rank(i, false), bv2.rank(i, false)) << V(i);
  }
  for(uint32_t i = 0; i < count[0]; i++) {
    ASSERT_EQ(bv.select(i + 1, false), bv2.select(i + 1, false)) << V(i); //select is 1-indexed
  }
  for(uint32_t i = 0; i < count[1]; i++) {
    ASSERT_EQ(bv.select(i + 1, true), bv2.select(i + 1, true)) << V(i);
  }
  // ASSERT_EQ(bv.select(count[0] + 2, false), bv.size());
  // ASSERT_EQ(bv.select(count[1] + 2, true), bv.size());
}

TEST(TestDynamicBV, TestRemoveBalanceLeftA) {
  // (128, 128) -> (64, 3 * 64 + 2) remove left, balance left a)
  BV<4> bv; NaiveBV bv2;
  uint32_t N = 64 + 3 * 64 + 2;
  for(uint32_t i = 0; i < 2 * 128; i++) {
    bool b = (i % 3) == 0;
    bv.insert(0, b); bv2.insert(0, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t i = 0; i < 64 + 2; i++) {
    bool b = (i % 3) == 0;
    bv.insert(150, b); bv2.insert(150, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 64; j++) {
    bv.remove(0); bv2.remove(0);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < N - (1 + j); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    } 
  }
}

TEST(TestDynamicBV, TestRemoveBalanceLeftB) {
  // (128 (128  128)) remove left, balance left b)
  BV<4> bv; NaiveBV bv2;
  uint32_t N = 3 * 128;
  for(uint32_t i = 0; i < N; i++) {
    bool b = (i % 2) == 0;
    bv.insert(i, b); bv2.insert(i, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  bv.insert(150, true); bv2.insert(150, true); //avoid merge

  for(uint32_t j = 0; j < 128; j++) {
    bv.remove(0); bv2.remove(0);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < N  + 1 - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveBalanceRightA) {
  // 120 | 64 -> remove right, balance right a)
  BV<2> bv; NaiveBV bv2;
  for(uint32_t i = 0; i < 120 + 64; i++) {
    bool b = (i % 2) == 0;
    bv.insert(0, b);
    bv2.insert(0, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 64; j++) {
    bv.remove(120);
    bv2.remove(120);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < 120 + 64 - (i + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveBalanceRightB) {
  // ((128, 250), 128) -> ((128, 250), 64) -> ((128, 192), 122), alignment is called, remove right, balance right b)
  BV<4> bv; NaiveBV bv2;
  uint32_t N = 128 + 250 + 128;
  for(uint32_t i = 0; i < 3 * 128; i++) {
    bool b = (i % 2) == 0;
    bv.insert(0, b); bv2.insert(0, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t i = 0; i < 250 - 128; i++) {
    bool b = (i % 2) == 0;
    bv.insert(130, b); bv2.insert(130, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 64; j++) {
    bv.remove(128 + 250 + 10); bv2.remove(128 + 250 + 10);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < N - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveMergeLeftARoot) {
  // 64 | 64 -> remove both sides, //merge left a)
  BV<2> bv; NaiveBV bv2;
  for(uint32_t i = 0; i < 128; i++) {
    bool b = (i % 2) == 0;
    bv.insert(0, b);
    bv2.insert(0, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 127; j++) {
    uint32_t k = (j % 2 == 0)? 0 : 127 - j;
    bv.remove(k);
    bv2.remove(k);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < 128 - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveMergeLeftAGrandparent) {
  // ((64 64) 64) -> shrink 1. & 2. block,   //merge left a)
  BV<2> bv; NaiveBV bv2;
  for(uint32_t i = 0; i < 3 * 64; i++) {
    bool b = (i % 5) == 0;
    bv.insert(0, b); bv2.insert(0, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 127; j++) {
    uint32_t k = (j % 2 == 0)? 0 : 127 - j;
    bv.remove(k); bv2.remove(k);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < 3 * 64 - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveMergeLeftB) {
  // (64 (64 64)) -> shrink 1. & 2. block,   //merge left b)
  BV<2> bv; NaiveBV bv2;
  for(uint32_t i = 0; i < 3 * 64; i++) {
    bool b = (i % 5) == 0;
    bv.insert(i, b); bv2.insert(i, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 127; j++) {
    uint32_t k = (j % 2 == 0)? 0 : 127 - j;
    bv.remove(k); bv2.remove(k);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < 3 * 64 - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveMergeRightA) {
  // ((64 64) 64) -> shrink 1. & 2. block,   //merge right a)
  BV<2> bv; NaiveBV bv2;
  for(uint32_t i = 0; i < 3 * 64; i++) {
    bool b = (i % 4) == 0;
    bv.insert(0, b); bv2.insert(0, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 127; j++) {
    uint32_t k = (j % 2 == 0)? 127 - j : 0;
    bv.remove(k); bv2.remove(k);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < 3 * 64 - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveMergeRightB) {
  // ((64 64) 64) -> ((64 48) 64) -> ((64 48) 16) (64 64) -> shrink 2. & 3. block,   //merge right b)
  BV<2> bv; NaiveBV bv2;
  for(uint32_t i = 0; i < 3 * 64; i++) {
    bool b = (i % 4) == 0;
    bv.insert(0, b); bv2.insert(0, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < 16; j++) {
    bv.remove(64); bv2.remove(64);
  }
  for(uint32_t j = 0; j < 48 + 5; j++) {
    bv.remove(64 + 48); bv2.remove(64 + 48);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < 64 + 48 + 64 - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemove1) {
  BV<8> bv; NaiveBV bv2;
  uint32_t N = 128 * 64;
  for(uint32_t i = 0; i < N; i++) {
    uint32_t j = 0; bool b = (i % 3) == 0;
    bv.insert(j, b); bv2.insert(j, b);
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t j = 0; j < N; j++) {
    bv.remove(0); bv2.remove(0);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < N - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}

TEST(TestDynamicBV, TestRemoveMid) {
  BV<8> bv; NaiveBV bv2;
  uint32_t N = 10000;
  std::vector<bool> seq = random_gen::random_bit_string(N);
  for(bool b : seq) {
    bv.insert(0, b); bv2.insert(0, b);
  }
  for(uint32_t j = 0; j < N; j++) {
    ASSERT_EQ(bv.size(), bv2.size());
    bv.remove(bv.size() / 2); 
    bv2.remove(bv2.size() / 2);
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t i = 0; i < N - (j + 1); i++) {
      ASSERT_EQ(bv.access(i), bv2.access(i)) << V(i);
    }
  }
}


TEST(TestDynamicBV, RandomOp) {
  BV<8> bv; NaiveBV bv2;
  //random is always same sequence
  uint32_t N = 10000;
  for(uint32_t i = 0; i < N; i++) {
    uint32_t j = random_gen::random_index(i + 1);
    bool b = random_gen::random_bool();
    bv.insert(j, b);
    bv2.insert(j, b);
    ASSERT_EQ(bv.size(), bv2.size());
    ASSERT_EQ(bv.check_invariant(), true);
  }
  for(uint32_t i = 0; i < N * 2; i++) {
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
    ASSERT_EQ(bv.check_invariant(), true);
  }

  for(uint32_t i = 0; i < N; i++) {
    ASSERT_EQ(bv.access(i), bv2.access(i));
  }
}

//disabled, to test with DynamicSmallBV<N>, which dont implement this method
TEST(TestDynamicBV, DISABLED_TestFastConstructionOffset) {
  uint32_t N = 2 * 6400 + 10;
  std::vector<bool> seq = random_gen::random_bit_string(N);
  std::vector<uint64_t> seq_compressed = utility::convert_bitstring(seq);
  NaiveBV bv2(seq);
  for(uint32_t off = 0; off < 10; off++) {
    uint32_t num_bits = N - off;
    BV<8> bv(seq_compressed, num_bits); 
    for(uint32_t j = 0; j < num_bits; j++) {
      ASSERT_EQ(bv.access(j), bv2.access(j));
      ASSERT_EQ(bv.rank(j, true), bv2.rank(j, true));
    }
  }
}

TEST(TestDynamicBV, TestFastConstructionSizes) {
  uint32_t N = 64 * 1000;
  std::vector<bool> seq = random_gen::random_bit_string(N);
  std::vector<bool> sub_seq;
  std::vector<uint32_t> sizes = {30, 1045, 3056, 8193, 12429, 55555, 63999};
  std::vector<uint64_t> seq_compressed;
  NaiveBV bv2(seq);
  for(auto num_bits : sizes) {
    ASSERT_EQ(num_bits < N, true);
    sub_seq.assign(seq.begin(), seq.begin() + num_bits);
    seq_compressed = utility::convert_bitstring(sub_seq);
    BV<2> bv(seq_compressed, num_bits); 
    ASSERT_EQ(bv.check_invariant(), true);
    for(uint32_t j = 0; j < num_bits; j++) {
      ASSERT_EQ(bv.access(j), bv2.access(j));
      ASSERT_EQ(bv.rank(j, true), bv2.rank(j, true));
    }
  }
}