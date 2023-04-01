#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>
#include <stack>
#include <utility>
#include <string>
#include <algorithm>
#include <cmath>
#include <variant>

#include "src/static_small_bv.h"

//debug
#define _ << " " <<
#define NL "\n"
#define V(x) std::string(#x "=") << (x) << " " //"x=...

/* Dynamic bitvector implementation with std::variant instead of tagged pointers. */

//N = number of uint64_t in leaf, N must be even so that split function works (can split between words)
template<unsigned int N = 2> 
class DynamicBVVariant {

    private:
        struct InnerNode;
        
    public:
        using int_type = uint32_t;
        using Leaf = StaticSmallBV<N>;
        using NodePtr = std::variant<InnerNode*, Leaf*>;

    private:
        struct InnerNode {
            NodePtr left;
            NodePtr right;
            int_type height{0};
            int_type num{0};
            int_type ones{0};
            ~InnerNode() {
                if(InnerNode** pt = std::get_if<InnerNode*>(&left)) delete *pt;
                else if(Leaf** pt = std::get_if<Leaf*>(&left)) delete *pt;
                if(InnerNode** pt = std::get_if<InnerNode*>(&right)) delete *pt;
                else if(Leaf** pt = std::get_if<Leaf*>(&right)) delete *pt;
            }

            InnerNode() {
                left = NULL; right = NULL; height = 1; num = 0; ones = 0;
            }

            //use after split
            InnerNode(Leaf* l, Leaf* r) {
                left = l; right = r; height = 1;
                num = l -> size(); ones = l -> rank(num, true);
            }
        };

    public:
        DynamicBVVariant() {
            root = new Leaf(); //start with only a leaf
        }

        ~DynamicBVVariant() {
            if(InnerNode** pt = std::get_if<InnerNode*>(&root)) delete *pt;
            else if(Leaf** pt = std::get_if<Leaf*>(&root)) delete *pt;
        }
        
        bool access(int_type i) {
            NodePtr cur = root;
            while(std::holds_alternative<InnerNode*>(cur)) {
                InnerNode* node = std::get<InnerNode*>(cur);
                if(i < node -> num) {cur = node -> left;}
                else {
                    i -= node -> num; 
                    cur = node -> right;
                }
            }
            return std::get<Leaf*>(cur) -> access(i);
        }

        void flip(int_type i) {
            NodePtr cur = root;
            std::stack<InnerNode*> stack;
            while(std::holds_alternative<InnerNode*>(cur)) {
                InnerNode* node = std::get<InnerNode*>(cur);
                if(i < node -> num) {
                    cur = node -> left; 
                    stack.push(node); //only push "left" nodes for update
                }
                else {
                    i -= node -> num;
                    cur = node -> right;
                }
            }
            Leaf* leaf = std::get<Leaf*>(cur);
            leaf -> flip(i);
            int32_t diff = (leaf -> access(i))? 1 : -1;
            while(!stack.empty()) {
                InnerNode* node = stack.top(); stack.pop();
                node -> ones += diff; 
            }
        }

        void insert(int_type i, bool b) {
            NodePtr cur = root;
            std::stack<std::pair<InnerNode*, bool>> stack;
            bool left_right = false;
            while(std::holds_alternative<InnerNode*>(cur)) {
                InnerNode* node = std::get<InnerNode*>(cur);
                if(i <= node -> num) {
                    node -> num++;
                    node -> ones += b;
                    cur = node -> left;
                    left_right = LEFT;
                }
                else {
                    i -= node -> num;
                    cur = node -> right;
                    left_right = RIGHT;
                }
                stack.push({node, left_right});
            }
            Leaf* leaf = std::get<Leaf*>(cur);
            leaf -> insert(i, b);
            if(leaf -> is_full()) {
                Leaf *r_leaf = new Leaf();
                leaf -> split(*r_leaf);
                InnerNode *splitter = new InnerNode(leaf, r_leaf);
                if(stack.empty()) { root = splitter;} //root pointed to leaf, now to splitter
                else {
                    auto[parent, dir] = stack.top();
                    if(dir == LEFT) parent -> left = splitter;
                    else            parent -> right = splitter;
                }
                InnerNode* child = NULL;
                while(!stack.empty()) { 
                    auto[node, path_dir] = stack.top();
                    if(child != NULL && is_imbalanced(child)) {
                        if(path_dir == LEFT) node -> left  = rebalance(child);
                        else                 node -> right = rebalance(child);
                    } 
                    update_height(node); //first rebalance than update height
                    child = node;
                    stack.pop(); 
                }
                if(child != NULL && is_imbalanced(child)) {
                    root = rebalance(child); 
                }
            }
        }

        int_type rank(int_type i, bool b) {
            NodePtr cur = root;
            int_type ones = 0;
            int_type inital_i = i;
            while(std::holds_alternative<InnerNode*>(cur)) {
                InnerNode* node = std::get<InnerNode*>(cur);
                if(i < node -> num) {cur = node -> left;}
                else {
                    i -= node -> num; ones += node -> ones;
                    cur = node -> right;
                }
            }
            ones +=  std::get<Leaf*>(cur) -> rank(i, true);
            return b * ones + (b^1) * (inital_i - ones);
        }

        int_type select(int_type i, bool b) {
            assert(i > 0); //beginn at 1.
            NodePtr cur = root;
            int_type pos = 0;
            while(std::holds_alternative<InnerNode*>(cur)) {
                InnerNode* node = std::get<InnerNode*>(cur);
                int_type zeros_or_ones = b? (node -> ones) : (node -> num) - (node -> ones);
                if(i <= zeros_or_ones) {cur = node -> left;} // <= , select is 1-indexed
                else {
                    i -= zeros_or_ones;
                    pos += node -> num;
                    cur = node -> right;
                }
            }
            Leaf* leaf = std::get<Leaf*>(cur);
            int_type local_select = leaf -> select(i, b);
            assert(local_select != leaf -> max_size()); //did not found index
            return pos + local_select;
        }

        bool check_invariant() {
            return check_invariant_rec(root);
        }

        void pr_tree() {
            pr_tree_rec(root, 0);
        }

    private:
        bool check_invariant_rec(NodePtr v) {
            if(std::holds_alternative<Leaf*>(v)) return true;
            InnerNode* node = std::get<InnerNode*>(v);
            int32_t balance = get_balance(node);
            if(std::abs(balance) > 1) return false;
            if(!check_invariant_rec(node -> left)) return false;
            if(!check_invariant_rec(node -> right)) return false;
            return true;
        }

        void update_height(InnerNode* v) { 
            v -> height = 1 + std::max(get_height(v -> left), get_height(v -> right));
        }

        int_type get_height(NodePtr v) {
            if(std::holds_alternative<Leaf*>(v)) return 0;
            return std::get<InnerNode*>(v) -> height;
        }

        int32_t get_balance(InnerNode* v) {
            return get_height(v -> left) -  get_height(v -> right);
        }

        bool is_imbalanced(InnerNode* v) {
            return std::abs(get_balance(v)) > 1;
        }

        /*
            y                               x
           / \     Right Rotation          /  \
          x   T3   - - - - - - - >        T1   y 
         / \       < - - - - - - -            / \
        T1  T2     Left Rotation            T2  T3
        */
        InnerNode* right_rotate(InnerNode* y) {
            InnerNode *x =  std::get<InnerNode*>(y -> left);
            NodePtr T2 = x -> right;

            //rotation
            x -> right = y;
            y -> left = T2;

            //update ones
            y -> ones -= x -> ones;
            y -> num -= x -> num;

            update_height(y);
            update_height(x);
            return x; //new root
        }

        InnerNode* left_rotate(InnerNode* x) {
            InnerNode* y = std::get<InnerNode*>(x -> right);
            NodePtr T2 = y -> left;
            
            //rotation
            x -> right = T2;
            y -> left = x; 

            //update ones
            y -> ones += x -> ones;
            y -> num += x -> num;

            update_height(x);
            update_height(y);
            return y; //new root
        }

        //returns new subtree root
        InnerNode* rebalance(InnerNode* z) { 
            assert(is_imbalanced(z)); //otherwise returns NULL
            int32_t balance = get_balance(z);
            InnerNode *subtree_root = NULL;
            if(balance > 1) { //left
                InnerNode* y = std::get<InnerNode*>(z -> left);
                int32_t balance_y = get_balance(y);
                if(balance_y > 0) { //left left
                    subtree_root = right_rotate(z);
                }
                else { //left right
                    z -> left = left_rotate(y);
                    subtree_root = right_rotate(z);
                }
                /* left left
                         z                                      y 
                        / \                                   /   \
                       y   T4      Right Rotate (z)          x      z
                      / \          - - - - - - - - ->      /  \    /  \ 
                     x   T3                               T1  T2  T3  T4
                    / \
                  T1   T2

                */

                /* left right
                     z                               z                           x
                    / \                            /   \                        /  \ 
                   y   T4  Left Rotate (y)        x    T4  Right Rotate(z)    y      z
                  / \      - - - - - - - - ->    /  \      - - - - - - - ->  / \    / \
                T1   x                          y    T3                    T1  T2 T3  T4
                    / \                        / \
                  T2   T3                    T1   T2

                */
                
            }
            else if(balance < -1){
                InnerNode* y = std::get<InnerNode*>(z -> right);
                int32_t balance_y = get_balance(y);
                if(balance_y < 0) { //right right
                    subtree_root = left_rotate(z);
                }
                else {
                    z -> right = right_rotate(y);
                    subtree_root = left_rotate(z);
                }
                /* right right
                     z                                y
                    /  \                            /   \ 
                   T1   y     Left Rotate(z)       z      x
                       /  \   - - - - - - - ->    / \    / \
                     T2   x                     T1  T2 T3  T4
                         / \
                       T3  T4
                */

                /* right left
                       z                            z                            x
                      / \                          / \                          /  \ 
                    T1   y   Right Rotate (y)    T1   x      Left Rotate(z)   z      y
                        / \  - - - - - - - - ->     /  \   - - - - - - - ->  / \    / \
                       x   T4                      T2   y                  T1  T2  T3  T4
                      / \                              /  \
                     T2   T3                           T3   T4
                */
            }
            return subtree_root;
        }

        void pr_tree_rec(NodePtr v, int_type level) {
            if(std::holds_alternative<Leaf*>(v)) {
                std::cout << std::string(level + 1, ' ' ) << "[" << std::get<Leaf*>(v) -> size() << "]\n";
                // std::get<Leaf*>(v) -> pr_all();
                return;
            }
            InnerNode* node = std::get<InnerNode*>(v);
            std::cout << std::string(level, ' ' ) << "( height, num, ones:" 
            << " " << node ->height 
            << " " << node ->num 
            << " " << node ->ones 
            << "\n"; 
            pr_tree_rec(node -> left, level + 1);
            pr_tree_rec(node -> right, level + 1);
            std::cout << std::string(level, ' ' ) << ")" << "\n";
        }

        NodePtr root;
        uint64_t cur_size = 0;
        const static bool LEFT = false;
        const static bool RIGHT = true;
};
