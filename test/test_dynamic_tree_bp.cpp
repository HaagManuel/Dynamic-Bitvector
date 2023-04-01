#include <gtest/gtest.h>
#include <vector>

#include "src/dynamic_bv.h"
#include "src/naive_bv.h"
#include "src/random_generator.h"

#include "src/naive_tree_bp.h"
#include "src/dynamic_tree_bp.h"

#define BV DynamicBV
// #define BV DynamicBVOop
// #define BV DynamicBVVariant


NaiveTreeBP insert_random_leafs(uint n) {
  NaiveTreeBP bp; // (( ))
  bp.insert(0,0);
  bp.insert(0,0);
  bp.insert(0,1);
  bp.insert(0,1);
  bp.print();
  for(uint i = 0; i < n; i++) {
    uint j = random_gen::random_index(bp.size());
    while(!bp.access(j)) j = random_gen::random_index(bp.size()); //~ 2x
    bp.insert_child(j, 0, 0);
    std::cout << V(bp.size()) << "\n";
    bp.print();
  }
  assert(bp.is_parentheses_expr());
  return bp;
}

//root with n childs of pattern pat
NaiveTreeBP instance1(uint n, std::vector<bool> &pat) {
  std::vector<bool> v;
  v.push_back(1);
  for(uint i = 0; i < n; i++) {
    v.insert(v.end(), pat.begin(), pat.end());
  }
  v.push_back(0);
  NaiveTreeBP bp(v);
  assert(bp.is_parentheses_expr());
  return bp;
}

TEST(TestDynamicBP, TreeFromLecture) {
  //   0123456789012345678901
  //   abbcddeffggechhijjkkia 
  //   (()(()(()()))()(()()))
  //   (()( ()(( )()) )()( ()() ))
  std::vector<bool> v = {1,1,0,1, 1,0,1,1, 0,1,0,0, 0,1,0,1, 1,0,1,0, 0,0}; 
  NaiveTreeBP bp1(v);
  DynamicTreeBP<2> bp2;
  for(uint i = 0; i < v.size(); i++) {
    bp2.insert(i, v[i]);
  }
  for(uint i = 0; i < v.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
    if(bp1.access(i)) { //is open
       ASSERT_EQ(bp1.close(i), bp2.close(i)) << i << "\n";
    }
    else {
      ASSERT_EQ(bp1.open(i), bp2.open(i)) << i << "\n";
    }
  }
}

TEST(TestDynamicBP, TestClose1) {
  int N = 400;
  std::vector<bool> pat = {1,1,1,0,1,0,0,0};
  NaiveTreeBP bp1 = instance1(N, pat);
  DynamicTreeBP<2> bp2;
  for(int i = 0; i < bp1.size(); i++) {
    bp2.insert(i, bp1.access(i));
  }
  for(int i = 0 + 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
    if(bp1.access(i)) { //is "("
       ASSERT_EQ(bp1.close(i), bp2.close(i)) << i << "\n";
    }
  }
}

TEST(TestDynamicBP, TestOpen1) {
  // int N = 20;
  int N = 400;
  std::vector<bool> pat = {1,1,1,0,1,0,0,0};
  NaiveTreeBP bp1 = instance1(N, pat);
  DynamicTreeBP<2> bp2;
  for(int i = 0; i < bp1.size(); i++) {
    bp2.insert(i, bp1.access(i));
  }
  // bp2.pr_tree();
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
    if(!bp1.access(i)) { //is ")"
      // std::cout << "open " << i << "\n";
      ASSERT_EQ(bp1.open(i), bp2.open(i)) << i << "\n";
    }
  }
}

TEST(TestDynamicBP, TestEnclose1) {
  // int N = 20;
  int N = 400;
  std::vector<bool> pat = {1,1,1,0,1,0,0,0};
  NaiveTreeBP bp1 = instance1(N, pat);
  DynamicTreeBP<2> bp2;
  for(int i = 0; i < bp1.size(); i++) {
    bp2.insert(i, bp1.access(i));
  }
  // bp2.pr_tree();
  for(int i = 1; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
    if(bp1.access(i)) { //is "("
      // std::cout << "enclose " << i << "\n";
      ASSERT_EQ(bp1.enclose(i), bp2.enclose(i)) << i << "\n";
    }
  }
}

TEST(TestDynamicBP, TestInsertLeafsRoot1) {
  int N = 1000;
  NaiveTreeBP bp1;
  DynamicTreeBP<2> bp2;
  bp1.only_root();
  bp2.only_root();
  for(int i = 0; i < N; i++) {
    bp1.insert_child(0, 1, 0); //insert leaf at first position of root, 0 is always node 0
    bp2.insert_child(0, 1, 0);
  }
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
  }
}

//insert at second position
TEST(TestDynamicBP, TestInsertLeafsRoot2) {
  int N = 1000;
  NaiveTreeBP bp1;
  DynamicTreeBP<2> bp2;
  bp1.only_root();
  bp2.only_root();
  for(int i = 0; i < N; i++) {
    bp1.insert_child(0, 2, 0);
    bp2.insert_child(0, 2, 0);
  }
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
  }
}

//insert at last position
TEST(TestDynamicBP, TestInsertLeafsRoot3) {
  int N = 500;
  NaiveTreeBP bp1;
  DynamicTreeBP<2> bp2;
  bp1.only_root();
  bp2.only_root();
  for(int i = 0; i < N; i++) {
    // std::cout << V(i) << "\n";
    // bp1.print();
    // bp2.pr_tree();
    bp1.insert_child(0, i + 1, 0);
    bp2.insert_child(0, i + 1, 0);
  }
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
  }
}

//insert with different k
TEST(TestDynamicBP, TestInsertLeafsRoot4) {
  int N = 200;
  for(int k = 0; k <= N; k++) {
    NaiveTreeBP bp1;
    DynamicTreeBP<2> bp2;
    bp1.only_root();
    bp2.only_root();
    for(int i = 0; i < N; i++) {
      bp1.insert_child(0, i + 1, 0);
      bp2.insert_child(0, i + 1, 0);
    }
    bp1.insert_child(0, 1, k);
    bp2.insert_child(0, 1, k);
    for(int i = 0; i < bp1.size(); i++) {
      ASSERT_EQ(bp1.access(i), bp2.access(i));
    }
  }
}


//not needed in 2. part 
// TEST(TestDynamicBP, TestRecConstruction)

TEST(TestDynamicBP, TestInsertRandomBasicOp) {
  int N = 1000;
  NaiveTreeBP bp1;
  DynamicTreeBP<4> bp2;
  bp1.only_root();
  bp2.only_root();
  for(int i = 0; i < N; i++) {
    int32_t k = 1 + random_gen::random_index(i + 1); //random node, node number start at 0
    int32_t v = bp2.select(k, true);
    ASSERT_EQ(bp1.select(k, true), v);
    int32_t degree = bp2.children(v);
    ASSERT_EQ(bp1.children(v), degree);
    int32_t i1 = 1 + random_gen::random_index(degree + 1); //children number start at 1
    int32_t i2 = random_gen::random_index(degree + 2 - i1);  // + 2?
    bp1.insert_child(v, i1, i2);
    bp2.insert_child(v, i1, i2);
  }
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
    if(bp1.access(i)) { // (
      ASSERT_EQ(bp1.close(i), bp2.close(i));
      ASSERT_EQ(bp1.subtree_size(i), bp2.subtree_size(i));
      if(i != 0) {
        ASSERT_EQ(bp1.parent(i), bp2.parent(i));
      }
      int32_t degree = bp2.children(i);
      ASSERT_EQ(bp1.children(i), degree);
      for(int j = 1; j <= degree; j++) {
        ASSERT_EQ(bp1.child_i(i, j), bp2.child_i(i, j));
      }
    }
    else { // )
      ASSERT_EQ(bp1.open(i), bp2.open(i));
    }
  }
  ASSERT_EQ(bp2.check_node_excess(), true);
  ASSERT_EQ(bp2.check_invariant(), true);
}

//big random test
TEST(TestDynamicBP, TestInsertRandom1) {
  // int N = 15000, M = N;
  int N = 1500, M = N;
  ASSERT_EQ(M <= N, true);
  NaiveTreeBP bp1;
  DynamicTreeBP<4> bp2; 
  bp1.only_root();
  bp2.only_root();
  for(int i = 0; i < N; i++) {
    // bp1.print();
    // bp2.pr_tree();
    int32_t k = 1 + random_gen::random_index(i + 1); //random node, node number start at 0
    int32_t v = bp2.select(k, true);
    ASSERT_EQ(bp1.select(k, true), v);
    int32_t degree = bp2.children(v);
    ASSERT_EQ(bp1.children(v), degree);
    int32_t i1 = 1 + random_gen::random_index(degree + 1); //children number start at 1
    int32_t i2 = random_gen::random_index(degree + 2 - i1);  // + 2?
    // std::cout << V(i) _ V(v) _ V(i1) _ V(i2) << "\n";
    bp1.insert_child(v, i1, i2);
    bp2.insert_child(v, i1, i2);
  }
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
  }
  ASSERT_EQ(bp2.check_node_excess(), true);

  // std::cout << "deleting nodes \n";
  for(int i = 0; i < M; i++) {
    // std::cout << V(i) << "\n";
    int32_t k = 1 + random_gen::random_index(N + 1 - i);
    int32_t v = bp2.select(k, true);
    ASSERT_EQ(bp1.select(k, true), v);
    ASSERT_EQ(bp1.close(v), bp2.close(v));
    bp1.delete_node(v);
    bp2.delete_node(v);
    // bp1.print();
    // bp2.pr_tree();
    ASSERT_EQ(bp2.check_node_excess(), true);
    for(int i = 0; i < bp1.size(); i++) {
      ASSERT_EQ(bp1.access(i), bp2.access(i)) << V(i) << V(bp1.size());
    }
  }
}


TEST(TestDynamicBP, TestDelete1) {
  // int N = 20000, M = 1000, childs = 5;
  // int N = 10000, childs = 3, M = N * childs;
  int N = 1000, childs = 3, M = N * childs;
  NaiveTreeBP bp1;
  DynamicTreeBP<4> bp2; 
  bp1.only_root();
  bp2.only_root();
  std::cout << "inserting nodes \n";
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < childs; j++) {
      // std::cout << V(i) _ V(j) << "\n";
      bp1.insert_child(i, 1, 0); //only leafs at first position
      bp2.insert_child(i, 1, 0);
      // ASSERT_EQ(bp2.check_invariant(), true);
      // ASSERT_EQ(bp2.check_node_excess(), true);
      // for(int i = 0; i < bp1.size(); i++) {
      //   ASSERT_EQ(bp1.access(i), bp2.access(i));
      // }
    }
  }
  ASSERT_EQ(bp2.check_invariant(), true);
  ASSERT_EQ(bp2.check_node_excess(), true);
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
  }

  std::cout << "deleting nodes \n";
  for(int i = 0; i < M; i++) {
    // std::cout << V(i) << "\n";
    ASSERT_EQ(bp1.size(), bp2.size());
    int32_t k = 1 + random_gen::random_index(bp1.size() / 2);
    int32_t v = bp1.select(k, true);
    ASSERT_EQ(bp2.select(k, true), v);
    ASSERT_EQ(bp1.access(v), true);
    
    // bp1.print();
    // bp2.pr_tree();
    // std::cout << V(bp1.size()) _ V(bp2.size()) _ V(v) << "\n";
    // std::cout << V(bp1.close(v)) << "\n";
    // std::cout << V(bp2.close(v)) << "\n";

    ASSERT_EQ(bp1.close(v), bp2.close(v));
    bp1.delete_node(v);
    bp2.delete_node(v);
    ASSERT_EQ(bp2.check_node_excess(), true);
    for(int i = 0; i < bp1.size(); i++) {
      ASSERT_EQ(bp1.access(i), bp2.access(i)) << V(i) << V(bp1.size());
    }
  }
}


TEST(TestDynamicBP, TestOutDegree) {
  int N = 5000;
  NaiveTreeBP bp1;
  DynamicTreeBP<4> bp2; 
  bp1.only_root();
  bp2.only_root();
  for(int i = 0; i < N; i++) {
    int32_t k = 1 + random_gen::random_index(i + 1); //random node, node number start at 0
    int32_t v = bp2.select(k, true);
    ASSERT_EQ(bp1.select(k, true), v);
    int32_t degree = bp2.children(v);
    ASSERT_EQ(bp1.children(v), degree);
    int32_t i1 = 1 + random_gen::random_index(degree + 1); //children number start at 1
    int32_t i2 = random_gen::random_index(degree + 2 - i1);  // + 2?
    // std::cout << V(i) _ V(v) _ V(i1) _ V(i2) << "\n";
    bp1.insert_child(v, i1, i2);
    bp2.insert_child(v, i1, i2);
  }
  for(int i = 0; i < bp1.size(); i++) {
    ASSERT_EQ(bp1.access(i), bp2.access(i));
  }
  ASSERT_EQ(bp2.check_node_excess(), true);
  
  std::vector<uint32_t> out_deg = bp2.all_out_deg();
  ASSERT_EQ(out_deg.size(), bp1.size() / 2);
  for(int32_t i = 1; i <= bp1.size() / 2; i++) {
      uint32_t v = bp1.select(i, true);
      uint32_t deg = bp2.children(v);
      ASSERT_EQ(bp1.children(v), deg);
      ASSERT_EQ(out_deg[i - 1], deg);
  }
}
