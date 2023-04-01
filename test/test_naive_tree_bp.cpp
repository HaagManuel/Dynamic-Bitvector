#include <gtest/gtest.h>

#include "src/naive_tree_bp.h"

TEST(TestNaiveBV, TreeFromLecture) {
  //   0123456789012345678901
  //   abbcddeffggechhijjkkia 
  //   (()(()(()()))()(()()))
  //   (()( ()(( )()) )()( ()() ))
  std::vector<bool> v = {1,1,0,1, 1,0,1,1, 0,1,0,0, 0,1,0,1, 1,0,1,0, 0,0}; 
  NaiveTreeBP bp(v);
  //closing
  EXPECT_EQ(bp.close(0), v.size() - 1); //a
  EXPECT_EQ(bp.close(1), 2); //b
  EXPECT_EQ(bp.close(3), 12); //c
  EXPECT_EQ(bp.close(4), 5); //d
  EXPECT_EQ(bp.close(6), 11); //e
  EXPECT_EQ(bp.close(7), 8); //f
  EXPECT_EQ(bp.close(9), 10); //g
  EXPECT_EQ(bp.close(13), 14); //h
  EXPECT_EQ(bp.close(15), 20); //i
  EXPECT_EQ(bp.close(16), 17); //j
  EXPECT_EQ(bp.close(18), 19); //k

  //open
  EXPECT_EQ(bp.open(v.size() - 1), 0); //a
  EXPECT_EQ(bp.open(2), 1); //b
  EXPECT_EQ(bp.open(12), 3); //c
  EXPECT_EQ(bp.open(5), 4); //d
  EXPECT_EQ(bp.open(11), 6); //e
  EXPECT_EQ(bp.open(8), 7); //f
  EXPECT_EQ(bp.open(10), 9); //g
  EXPECT_EQ(bp.open(14), 13); //h
  EXPECT_EQ(bp.open(20), 15); //i
  EXPECT_EQ(bp.open(17), 16); //j
  EXPECT_EQ(bp.open(19), 18); //k

  //enclose
  EXPECT_EQ(bp.enclose(1), 0); //b
  EXPECT_EQ(bp.enclose(3), 0); //c
  EXPECT_EQ(bp.enclose(4), 3); //d
  EXPECT_EQ(bp.enclose(6), 3); //e
  EXPECT_EQ(bp.enclose(7), 6); //f
  EXPECT_EQ(bp.enclose(9), 6); //g
  EXPECT_EQ(bp.enclose(13), 0); //h
  EXPECT_EQ(bp.enclose(15), 0); //i
  EXPECT_EQ(bp.enclose(16), 15); //j
  EXPECT_EQ(bp.enclose(18), 15); //k

  //childs
  EXPECT_EQ(bp.child_i(0,1), 1); //a -> b
  EXPECT_EQ(bp.child_i(0,2), 3); //a -> c
  EXPECT_EQ(bp.child_i(0,3), 13); //a -> h
  EXPECT_EQ(bp.child_i(0,4), 15); //a -> i
  EXPECT_EQ(bp.child_i(3,1), 4); //c -> d
  EXPECT_EQ(bp.child_i(3,2), 6); //c -> e
  EXPECT_EQ(bp.child_i(6,1), 7); //e -> f
  EXPECT_EQ(bp.child_i(6,2), 9); //e -> g
  EXPECT_EQ(bp.child_i(15,1), 16); //i -> j
  EXPECT_EQ(bp.child_i(15,2), 18); //i -> k

  //children
  EXPECT_EQ(bp.children(0), 4);
  EXPECT_EQ(bp.children(1), 0);
  EXPECT_EQ(bp.children(3), 2);
  EXPECT_EQ(bp.children(6), 2);
  EXPECT_EQ(bp.children(15), 2);
  EXPECT_EQ(bp.children(15), 2);
}
