#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>
#include <stack>
#include <utility>
#include <string>
#include <algorithm>
#include <cmath>
#include <tuple>

#include "src/static_small_bv.h"
#include "src/tagged_pointer.h"


#include "src/dynamic_small_bv.h"


//debug
#define _ << " " <<
#define NL "\n"
#define V(x) std::string(#x "=") << (x) << " " //"x=...

// #define DEBUG_OUT false

#define DBG if constexpr(DEBUG_OUT) std::cout

/* Main Dynamic bitvector implementation with AVL Tree.
   N = number of uint64_t in leaf, N must be even so that split function works (can split between words)
*/
template<unsigned int N = 2> 
class DynamicBV {
    public:
        using int_type = uint32_t;
        using height_type = uint32_t;
        using Leaf = StaticSmallBV<N>; //fixed size for leaf -> better runtime
        // using Leaf = DynamicSmallBV<N>; //reallocation leaf when size changes -> more space efficient, buggy at the moment

    private:
    static const bool DEBUG_OUT = false; //debug output
        struct Node {
            TaggedPointer left; // tag == true -> node, tag == false -> leaf
            TaggedPointer right;
            height_type height{0}; 
            int_type num{0};
            int_type ones{0};
            ~Node() {
                auto [leaf_pt_l, node_pt_l, mark_l] = cast_pointer(left);
                auto [leaf_pt_r, node_pt_r, mark_r] = cast_pointer(right);
                if(mark_l == IS_LEAF) { delete leaf_pt_l; }
                else {delete node_pt_l;}
                if(mark_r == IS_LEAF) { delete leaf_pt_r; }
                else {delete node_pt_r;}
            }

            Node() {
                left = NULL;
                right = NULL;
                height = 1;
                num = 0;
                ones = 0;
            }

            Node(TaggedPointer l, TaggedPointer r, int_type _num, int_type _ones) : left(l), right(r), num(_num), ones(_ones) {}

            //use after split
            Node(Leaf *l, Leaf *r) {
                left = TaggedPointer(l, IS_LEAF);
                right = TaggedPointer(r, IS_LEAF);
                height = 1;
                reinit_from_leaf(l);
            }

            void reinit_from_leaf(Leaf* l) {
                num = l -> size();
                ones = l -> num_ones();
            } 

            void invalidate_pointers() {
                left = NULL; right = NULL;
            }
        };

    public:
        DynamicBV() {
            //start with only a leaf
            Leaf *l = new Leaf();
            root = TaggedPointer(l, IS_LEAF);
        }

        ~DynamicBV() {
            auto [leaf_pt, node_pt, mark] = cast_pointer(root);
            if(mark == IS_LEAF) { delete leaf_pt; }
            else                { delete node_pt; }
        }
        
        DynamicBV(std::vector<uint64_t> &bits, uint32_t num_bits) {
            int32_t trailing_bits = 64 * bits.size() - num_bits; 
            assert(trailing_bits >= 0 && trailing_bits < 64);
            auto[root_pt, num, ones] = rec_construction(bits, 0, bits.size()); 
            root = root_pt;
            assert(check_invariant());
            //remove trailing bits
            TaggedPointer cur = root;
            while(is_node(cur)) cur = get_node(cur) -> right;
            get_leaf(cur) -> add_size(-trailing_bits); //decrease size of rightmost leaf
            //no updates of num and ones necessary since it is the rightmost leaf
            cur_size = num_bits;
        }

        //pointer from subtree, nums, ones,    start inclusive, end exclusive
        std::tuple<TaggedPointer, int_type, int_type> rec_construction(std::vector<uint64_t> &bits, uint32_t start, uint32_t end) {
            uint32_t len = end - start;
            uint32_t mid = (start + end) / 2;
            auto copy_to_leaf = [&](uint32_t from, uint32_t to) {
                Leaf* leaf = new Leaf();
                for(uint32_t i = from; i < to; i++) leaf -> set_word(bits[i], i - from);
                leaf -> set_size((to - from) * leaf -> BITS);
                return leaf;
            };
            if(len == N) {
                Leaf* l = copy_to_leaf(start, mid); Leaf* r = copy_to_leaf(mid, end);
                int_type l_num = l -> size();
                int_type l_ones = l -> num_ones();
                int_type r_num = r -> size();
                int_type r_ones = r -> num_ones();
                Node *node = new Node(l, r, l_num, l_ones);
                update_height(node);
                return {TaggedPointer(node, IS_NODE), l_num + r_num, l_ones + r_ones};
            }
            else if(len < N) {
                Leaf* leaf = copy_to_leaf(start, end);
                return {TaggedPointer(leaf, IS_LEAF), leaf -> size(), leaf -> num_ones()};
            }
            auto [l_pt, l_num, l_ones] = rec_construction(bits, start, mid);
            auto [r_pt, r_num, r_ones] = rec_construction(bits, mid, end);
            Node* node = new Node(l_pt, r_pt, l_num, l_ones);
            update_height(node);
            return {TaggedPointer(node, IS_NODE), l_num + r_num, l_ones + r_ones};
        }
        
        /*memory computation */
        struct MemoryData {
            uint64_t num_leafs;
            uint64_t num_nodes;
            uint64_t memory_leafs; //in bytes
            uint64_t memory_nodes; //in bytes
            uint64_t memory_total; //in bytes
            uint64_t memory_stored_bits; //in bytes
            uint64_t memory_all_bits; //in bytes
            double fraction_stored_used_bits; 
            double bits_per_stored_bit; 
            MemoryData(uint64_t _num_leafs, uint64_t _num_nodes, uint64_t _bits_set) {
                num_leafs = _num_leafs;
                num_nodes = _num_nodes;
                memory_leafs = sizeof(Leaf) * num_leafs;
                memory_nodes = sizeof(Node) * num_nodes;
                memory_total = memory_leafs + memory_nodes + sizeof(DynamicBV);
                memory_stored_bits = _bits_set / 8; //bit -> byte
                memory_all_bits = 8 * num_leafs * N; //8 byte per uint64_t per leaf
                fraction_stored_used_bits = (double) memory_stored_bits / memory_all_bits;
                bits_per_stored_bit = (double) memory_total / memory_stored_bits;
            }

            void print() {
                std::string header = std::string(50, '-');
                std::cout << header << NL
                          << "memory details, units in Bytes" << NL
                          << V(num_leafs) << NL
                          << V(num_nodes) << NL
                          << V(memory_leafs) << NL
                          << V(memory_nodes) << NL 
                          << V(memory_total) << NL
                          << V(memory_stored_bits) << NL
                          << V(memory_all_bits) << NL
                          << V(fraction_stored_used_bits) << NL
                          << V(bits_per_stored_bit) << NL
                          << header << NL
                        ;
            }
        };
        //inner_nodes, leafs, raw bits
        MemoryData memory_consumption() {
            auto[num_leafs, num_nodes] = leafs_and_nodes();
            return MemoryData(num_leafs, num_nodes, cur_size);
        }

        std::tuple<uint64_t, uint64_t> leafs_and_nodes() {
            return leafs_and_nodes_rec(root);
        }

        std::tuple<uint64_t, uint64_t> leafs_and_nodes_rec(TaggedPointer pt) {
            if(is_leaf(pt)) { return {1,0}; }
            else {
                Node* node = get_node(pt);
                auto[l1, n1] = leafs_and_nodes_rec(node -> left);
                auto[l2, n2] = leafs_and_nodes_rec(node -> right);
                return {l1 + l2, n1 + n2 + 1};
            }
        }
        /*memory computation */


        bool access(int_type i) {
            assert(i < cur_size);
            TaggedPointer cur = root;
            while(is_node(cur)) {
                Node* node = get_node(cur);
                if(i < node -> num) {cur = node -> left;}
                else {
                    i -= node -> num;
                    cur = node -> right;
                }
            }
            return get_leaf(cur) -> access(i);
        }


        /* flip variants */

        /* iterative optimized */
        void flip1(int_type i) {    
            assert(i < cur_size);
            TaggedPointer cur = root;
            std::stack<Node*> stack;
            while(is_node(cur)) {
                Node* node = get_node(cur);
                if(i < node -> num) {
                    cur = node -> left; 
                    stack.push(node); // only push "left" nodes
                }
                else {
                    i -= node -> num;
                    cur = node -> right;
                }
            }
            Leaf* leaf = get_leaf(cur);
            int32_t diff = (leaf -> access(i)) == 0? 1 : -1;
            leaf -> flip(i);
            while(!stack.empty()) {
                auto node = stack.top(); stack.pop(); //only update "left" nodes
                node -> ones += diff; 
            }
            return;
        }


        /* iterative */
        void flip2(int_type i) {
            assert(i < cur_size);
            TaggedPointer cur = root;
            std::stack<std::pair<Node*,bool>> stack;
            while(is_node(cur)) {
                Node* node = get_node(cur);
                if(i < node -> num) {
                    cur = node -> left; 
                    stack.push({node, LEFT});
                }
                else {
                    i -= node -> num;
                    cur = node -> right;
                    stack.push({node, RIGHT});
                }
            }
            Leaf* leaf = get_leaf(cur);
            int32_t diff = (leaf -> access(i)) == 0? 1 : -1;
            leaf -> flip(i);
            while(!stack.empty()) {
                auto [node,lr] = stack.top(); stack.pop();
                if(lr == LEFT) node -> ones += diff; 
            }
            return;
        }

        /* recursive -> best*/
        void flip(int_type i) {
            flip_rec(root, i);
        }

        bool flip_rec(TaggedPointer pt, int_type i)  { 
            auto [leaf, node, mark] = cast_pointer(pt);
            if(mark == IS_LEAF) {
                leaf -> flip(i);
                return leaf -> access(i);
            }
            bool b;
            if(i < node -> num) {
                b = flip_rec(node -> left, i);
                node -> ones += b? 1 : -1;
            }
            else { b = flip_rec(node -> right, i - node -> num);}
            return b;
        }

        /* flip variants */

        int_type rank(int_type i, bool b) {
            assert(i < cur_size + 1);
            TaggedPointer cur = root;
            int_type ones = 0;
            int_type inital_i = i;
            while(is_node(cur)) {
                Node *node = get_node(cur);
                if(i < node -> num) {cur = node -> left;}
                else {
                    i -= node -> num;
                    ones += node -> ones;
                    cur = node -> right;
                }
            }
            ones += get_leaf(cur) -> rank(i, true);
            return b * ones + (b^1) * (inital_i - ones);
        }

        int_type select(int_type i, bool b) {
            assert(i > 0); //beginn at 1.
            TaggedPointer cur = root;
            int_type pos = 0;
            while(is_node(cur)) {
                Node *node = get_node(cur);
                int_type zeros_or_ones = b? (node -> ones) : (node -> num) - (node -> ones);
                if(i <= zeros_or_ones) {cur = node -> left;} // <= , select is 1-indexed
                else {
                    i -= zeros_or_ones;
                    pos += node -> num;
                    cur = node -> right;
                }
            }
            Leaf* leaf = get_leaf(cur);
            int_type local_select = leaf -> select(i, b);
            if(local_select == leaf -> max_size()) {
                return cur_size; //did not found index
            }
            return pos + local_select;
        }

        void update_path(std::stack<std::pair<Node*, bool>> &stack) {
            Node *child = NULL;
            while(!stack.empty()) { 
                auto[node, path_dir] = stack.top();
                if(child != NULL && is_imbalanced(child)) {
                    TaggedPointer pt = TaggedPointer(rebalance(child), IS_NODE);
                    assign_pointer(pt, node, path_dir); //assign subtree from rebalancing
                } 
                //first rebalance than update height
                node -> height = 1 + std::max(get_height(node -> left), get_height(node -> right));
                child = node;
                stack.pop(); 
            }
            if(child != NULL && is_imbalanced(child)) {
                root = TaggedPointer(rebalance(child), IS_NODE); 
            }
        }

        void insert(int_type i, bool b) {
            assert(i < cur_size + 1);
            cur_size++;
            TaggedPointer cur = root;
            std::stack<std::pair<Node*, bool>> stack;
            bool left_right = false;
            while(is_node(cur)) {
                Node *node = get_node(cur);
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
            Leaf *leaf = get_leaf(cur);
            leaf -> insert(i, b);
            if(leaf -> is_full()) { //split if leaf is full
                TaggedPointer pt = split_leaf(leaf);
                if(stack.size() == 0) { root = pt;} //root pointed to leaf, now to splitter
                else {
                    auto[parent, lr] = stack.top();
                    assign_pointer(pt, parent, lr);
                }
                update_path(stack);
            }
        }

        TaggedPointer split_leaf(Leaf* leaf) {
            assert(leaf -> is_full());
            Leaf *r_leaf = new Leaf();
            leaf -> split(*r_leaf);
            Node *splitter = new Node(leaf, r_leaf);
            return TaggedPointer(splitter, IS_NODE);
        }

        void remove(int_type i) {
            assert(i < cur_size && cur_size > 0);
            cur_size--;
            TaggedPointer cur = root;
            std::stack<Node*> update_stack;
            std::stack<std::pair<Node*, bool>> stack;
            bool left_right = false;
            while(is_node(cur)) {
                Node *node = get_node(cur);
                if(i < node -> num) {
                    cur = node -> left;
                    left_right = LEFT;
                    update_stack.push(node);
                }
                else {
                    i -= node -> num;
                    cur = node -> right;
                    left_right = RIGHT;
                }
                stack.push({node, left_right});
            }
            Leaf *leaf = get_leaf(cur);
            bool removed_bit = leaf -> access(i);
            leaf -> remove(i);
    
            auto update_num_ones = [&]() {
                while(!update_stack.empty()) {
                    Node* node = update_stack.top(); update_stack.pop();
                    node -> num -= 1;
                    node -> ones -= removed_bit;
                }
            };
            if(is_leaf(root)) {  //do nothing
                update_num_ones();
                return;
            } 
            //root is splitter
            assert(stack.size() > 0); //splitter node
            Node* parent = NULL;
            Node* grandparent = NULL;
            bool lr_parent = LEFT, lr_grandparent = LEFT;
            std::tie(parent, lr_parent) = stack.top();
            stack.pop();
            if(stack.size() > 0) std::tie(grandparent, lr_grandparent) = stack.top();
            stack.push({parent, lr_parent}); 
            int_type MAX_SIZE = leaf -> MAX_SIZE;

            //first try to balance, than merge
            auto try_balance = [&](Leaf* a, Leaf* b) {
                Leaf::simple_balance(*a, *b);
                return true;
            };
            auto sum_small =  [&](Leaf* a, Leaf* b) { return a -> size() + b -> size() <= (MAX_SIZE / 2); };
            auto merge_leaf = [&](Leaf* target, Leaf* other) {
                target -> simple_append(*other);
            };
            bool merged = false;
            if(lr_parent == LEFT) {
                /* search right nghr
                 a)   parent             or     b)   parent
                     /    \                         /     \
                    leaf   leaf2                   leaf    inner1
                                                          /      \
                                                        leaf2   leaf3
                */
                TaggedPointer pt = parent -> right;
                if(is_leaf(pt)) { //a)
                    Leaf* leaf2 = get_leaf(pt);
                    if(sum_small(leaf, leaf2)) { 
                        DBG << "merge left a) \n";
                        merge_leaf(leaf, leaf2);
                        update_stack.pop();
                        stack.pop(); //delete parent from stack
                        TaggedPointer pt_leaf = TaggedPointer(leaf, IS_LEAF);
                        if(grandparent != NULL) assign_pointer(pt_leaf, grandparent, lr_grandparent);
                        else                    root = pt_leaf;
                        parent -> invalidate_pointers();
                        delete parent;
                        delete leaf2;
                        merged = true;
                    }
                    else if(try_balance(leaf, leaf2)) {
                        DBG << "balance left a) \n";
                        update_stack.pop(); //delete parent from update_stack, num and ones got updated already
                        parent -> reinit_from_leaf(leaf);
                    } 
                }
                else { // b)
                    Node* inner1 = get_node(pt);
                    Leaf* leaf2 = get_leaf(inner1 -> left);
                    if(sum_small(leaf, leaf2)) {
                        DBG << "merge left b) \n";
                        // leaf -> append_bv(*leaf2); //merge
                        merge_leaf(leaf, leaf2);
                        parent -> right = TaggedPointer(get_leaf(inner1 -> right), IS_LEAF);
                        update_stack.pop(); //delete parent from update_stack, num and ones got updated already
                        parent -> reinit_from_leaf(leaf);
                        inner1 -> invalidate_pointers();
                        delete inner1;
                        delete leaf2;
                        merged = true;
                    }
                    else if(try_balance(leaf, leaf2)) {
                        DBG << "balance left b) \n";
                        update_stack.pop(); //delete parent from update_stack, num and ones got updated already
                        parent -> reinit_from_leaf(leaf);
                        inner1 -> reinit_from_leaf(leaf2);
                    }  
                }
            }
            else  { //lr_parent == RIGHT
                /* search right nghr
                 a)   parent             or      b)  parent
                     /    \                         /     \
                    leaf2   leaf                  inner1    leaf
                                                 /      \
                                                leaf3   leaf2
                */
                TaggedPointer pt = parent -> left;
                if(is_leaf(pt)) { //a)
                    Leaf* leaf2 = get_leaf(pt);
                    if(sum_small(leaf2, leaf)) { //cannot force alignment!
                        DBG << "merge right a) \n";
                        merge_leaf(leaf2, leaf);
                        // update_stack.pop(); //parent is not on update stack, since it comes from right
                        stack.pop(); //delete parent from stack
                        TaggedPointer pt_leaf = TaggedPointer(leaf2, IS_LEAF);
                        if(grandparent != NULL) assign_pointer(pt_leaf, grandparent, lr_grandparent);
                        else                    root = pt_leaf;
                        parent -> invalidate_pointers();
                        delete parent;
                        delete leaf;
                        merged = true;
                    }
                    else if(try_balance(leaf2, leaf)) {
                        DBG << "balance right a) \n";
                        parent -> reinit_from_leaf(leaf2);
                    }
                }
                else { // b)
                    Node* inner1 = get_node(pt);
                    Leaf* leaf2 = get_leaf(inner1 -> right);
                    Leaf* leaf3 = get_leaf(inner1 -> left);
                    if(sum_small(leaf2, leaf)) { //cannot force alignment!
                        DBG << "merge right b) \n";
                        // leaf2 -> append_bv(*leaf); //merge
                        merge_leaf(leaf2, leaf);
                        parent -> left = TaggedPointer(leaf3, IS_LEAF);
                        parent -> right = TaggedPointer(leaf2, IS_LEAF);
                        parent -> reinit_from_leaf(leaf3);
                        inner1 -> invalidate_pointers();
                        delete inner1;
                        delete leaf;
                        merged = true;
                    } 
                    else if(try_balance(leaf2, leaf)) {
                        DBG << "balance right b) \n";
                        parent -> num = leaf3 -> size() + leaf2 -> size();
                        parent -> ones = leaf3 -> num_ones() + leaf2 -> num_ones();
                    }
                }
            }
            update_num_ones();
            if(merged) {
                update_path(stack); //rebalancing, update height
            }
        }

        //returns new subtree root
        Node* rebalance(Node* z) { 
            assert(is_imbalanced(z)); //otherwise returns NULL
            int32_t balance = get_balance(z);
            Node *subtree_root = NULL;
            if(balance > 1) { //left
                auto [leaf, node, mark] = cast_pointer(z -> left);
                assert(mark == IS_NODE);
                Node* y = node;
                int32_t balance_y = get_balance(y);
                if(balance_y > 0) { //left left
                    // std::cout << "left left \n";
                    assert(y -> left.getMark() == IS_NODE);
                    subtree_root = right_rotate(z);
                }
                else { //left right
                    // std::cout << "left right\n";
                    assert(y -> right.getMark() == IS_NODE);
                    z -> left = TaggedPointer(left_rotate(y), IS_NODE); 
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
                auto [leaf, node, mark] = cast_pointer(z -> right);
                assert(mark == IS_NODE);
                Node* y = node;
                int32_t balance_y = get_balance(y);
                if(balance_y < 0) { //right right
                    // std::cout << "right right \n";
                    subtree_root = left_rotate(z);
                }
                else {
                    // std::cout << "right left \n";
                    z -> right = TaggedPointer(right_rotate(y), IS_NODE);
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

        uint64_t size() { return cur_size; }

        uint64_t leaf_blocks() { return N;}

        void pr_tree() {
            pr_tree_rec(root, 0);
        }

        //|balance| <= 1
        bool check_invariant() {
            return check_invariant_rec(root);
        }

        uint32_t blocks_per_leaf() {return N;}

    private:
        void update_height(Node* v) { 
            v -> height = 1 + std::max(get_height(v -> left), get_height(v -> right));
        }

        bool check_invariant_rec(TaggedPointer pt) {
            auto[leaf, node, mark] = cast_pointer(pt);
            if(mark == IS_LEAF) return true;
            int32_t balance = get_balance(node);
            if(std::abs(balance) > 1) return false;
            if(!check_invariant_rec(node -> left)) return false;
            if(!check_invariant_rec(node -> right)) return false;
            return true;
        }
        /*
            y                               x
           / \     Right Rotation          /  \
          x   T3   - - - - - - - >        T1   y 
         / \       < - - - - - - -            / \
        T1  T2     Left Rotation            T2  T3
        */
        Node* right_rotate(Node* y) {
            auto [leaf, node, mark] = cast_pointer(y -> left);
            assert(mark == IS_NODE);
            Node *x = node;
            TaggedPointer t2 = x -> right;
            //rotation
            x -> right = TaggedPointer(y, IS_NODE);
            y -> left = t2;

            //update ones
            y -> ones -= x -> ones;
            y -> num -= x -> num;

            update_height(y);
            update_height(x);
            return x; //new root
        }

        Node* left_rotate(Node* x) {
            auto [leaf, node, mark] = cast_pointer(x -> right);
            assert(mark == IS_NODE);
            Node* y = node;
            TaggedPointer t2 = y -> left;
            
            //rotation
            x -> right = t2;
            y -> left = TaggedPointer(x, IS_NODE); 

            //update ones
            y -> ones += x -> ones;
            y -> num += x -> num;

            update_height(x);
            update_height(y);
            return y; //new root
        }

        void pr_tree_rec(TaggedPointer pt, int_type level) {
            auto [leaf_pt, node_pt, mark] = cast_pointer(pt);
            if(mark == IS_LEAF) {
                std::cout << std::string(level + 1, ' ' ) << "[" << leaf_pt -> size() << "]\n";
                // leaf_pt -> pr_all();
                return;
            }
            std::cout << std::string(level, ' ' ) << "( height, num, ones:" 
            << " " << node_pt ->height 
            << " " << node_pt ->num 
            << " " << node_pt ->ones 
            << "\n"; 
            pr_tree_rec(node_pt -> left, level + 1);
            pr_tree_rec(node_pt -> right, level + 1);
            std::cout << std::string(level, ' ' ) << ")" << "\n";
        }

        static std::tuple<Leaf*, Node*, bool> cast_pointer(TaggedPointer &pt) { 
            return {(Leaf*)pt.getRef(), (Node*)pt.getRef(), pt.getMark()};
        }

        static bool is_leaf(TaggedPointer &pt) {
            return pt.getMark() == IS_LEAF;
        }

        static bool is_node(TaggedPointer &pt) {
            return pt.getMark() == IS_NODE;
        }

        static Leaf* get_leaf(TaggedPointer &pt) {
            return (Leaf*)(pt.getRef());
        }

        static Node* get_node(TaggedPointer &pt) {
            return (Node*)(pt.getRef());
        }

        static void assign_pointer(TaggedPointer &pt, Node* v, bool lr) {
            if(lr == LEFT) v -> left = pt;
            else           v -> right = pt;
        }

        height_type get_height(TaggedPointer &pt) {
            auto[leaf, node, mark] = cast_pointer(pt);
            if(mark == IS_LEAF) {return 0;}
            return node -> height;
        }

        int32_t get_balance(Node* v) {
            return get_height(v -> left) -  get_height(v -> right);
        }

        bool is_imbalanced(Node* v) {
            return std::abs(get_balance(v)) > 1;
        }

        void pr_bit(uint64_t x) {
            std::bitset<64> bs(x);
            std::cout << bs << "\n";
        }

        TaggedPointer root;
        uint64_t cur_size = 0;
        const static bool IS_LEAF = false;
        const static bool IS_NODE = true;
        const static bool LEFT = false;
        const static bool RIGHT = true;
};
