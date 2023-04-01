#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>
#include <stack>
#include <utility>
#include <string>
#include <algorithm>
#include <cmath>

#include "src/oop/node.h"
#include "src/oop/inner_node.h"
#include "src/oop/leaf.h"

/* Dynamic bitvector implementation with inheritance for comparion run time overhead to tagged pointer variant.
   Remove operation is not implemented.
*/
template<unsigned int N = 2> 
class DynamicBVOop {
    public:
        using int_type = uint32_t;
        DynamicBVOop() {
            root = new Leaf<N>();
        }

        ~DynamicBVOop() {
            delete root;
        }


        bool access(int_type i) {return root -> access(i);}
        void flip(int_type i) {root -> flip(i);}
        void insert(int_type i , bool b) { root = root -> insert(i, b); cur_size++;}
        void remove(int_type) { /*not yet implemented*/ cur_size--;}
        
        int_type rank(int_type i , bool b) { return b? root -> rank1(i, 0) : i - root -> rank1(i, 0);}
        int_type select(int_type i , bool b) {return b? root -> select1(i, 0) : root -> select0(i, 0);}

        uint64_t size() {return cur_size;}
        uint64_t leaf_blocks() {return N;}

        bool check_invariant() {return root -> check_invariant();}
        void pr_tree() {
            root -> pr_tree(0);
        }
    private:
        Node* root;
        uint64_t cur_size{0};

};