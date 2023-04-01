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
#include <limits>

#include "src/static_small_bv.h"
#include "src/tagged_pointer.h"

//debug
#define _ << " " <<
#define NL "\n"
#define V(x) std::string(#x "=") << (x) << " " //"x=...


#define DBG if constexpr(DEBUG_OUT) std::cout

/* Main Dynamic Tree Balanced Parentheses implementation with AVL Tree.
   Code is copied and adapted from dynamic_bv.h.
   N = number of uint64_t in leaf, N must be even so that split function works (can split between words)
*/
template<unsigned int N = 2> 
class DynamicTreeBP {
    public:
        using int_type = uint32_t;
        // using height_type = uint8_t; //max heigthx 2^256 (requires correct balancing)
        using height_type = uint32_t; //no difference because of alignment
        using Leaf = StaticSmallBV<N>;

    private:
    static const bool DEBUG_OUT = false; //debug output

        struct Excess {
            Excess() {};

            Excess(Excess &l, Excess &r) {
                combine(l, r);
            }

            //used when scanning bits from left to right
            void update(bool bit) {
                e += bit? 1 : -1;
                if     (e < m)  m = e;
            }

            void update_right(bool bit) {
                e += bit? -1 : 1;
                if     (e < m_r)  m_r = e;
            }

            void combine(Excess &l, Excess &r) {
                e = l.e + r.e;
                m = std::min(l.m, l.e + r.m);
                m_r = std::min(l.m_r - r.e, r.m_r);
            }

            void print() {
                std::cout << V(e) _ V(m) _ V(m_r) << "\n";
            }

            int32_t min_excess_right_left() {
                return m - e;
            }

            int32_t e{0}; //excess left to right
            int32_t m{2}; //min local excess from left to right
            int32_t m_r{2}; //min local excess from right to left
        };

        struct Node {
            TaggedPointer left; // tag == true -> node, tag == false -> leaf
            TaggedPointer right;
            height_type height{0}; // -> maximum
            int_type num{0};
            int_type ones{0};
            Excess excess{};
            
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
                // recompute_excess();
            }

            void reinit_from_leaf(Leaf* l) {
                num = l -> size();
                ones = l -> num_ones();
            } 

            void invalidate_pointers() {
                left = NULL; right = NULL;
            }

            Excess excess_from_right() {
                Excess node_ex;
                node_ex.e = -excess.e; //excess left to right -> excess right to left
                node_ex.m = excess.m;
                node_ex.m_r = excess.m_r; 
                return node_ex;
            }
        };

    Excess excess_leaf(Leaf &l) {
        Excess ex;
        for(int_type i = 0; i < l.size(); i++) {
            ex.update(l.access(i));
        }
        return ex;
    }

    Excess excess_leaf_right(Leaf &l) {
        Excess ex;
        for(int32_t i = (int32_t)(l.size()) - 1; i >= 0; i--) {
            ex.update_right(l.access(i));
        }
        return ex;
    }

    Excess excess_node(TaggedPointer pt) {
        auto [leaf_pt, node_pt, mark] = cast_pointer(pt);
        if(mark == IS_LEAF) return excess_leaf(*leaf_pt);
        else                return node_pt -> excess;
    }

    Excess excess_node_right(TaggedPointer pt) {
        auto [leaf_pt, node_pt, mark] = cast_pointer(pt);
        if(mark == IS_LEAF) return excess_leaf_right(*leaf_pt);
        else                return node_pt -> excess_from_right();
    }

    void recompute_excess(Node *node) {
        Excess ex_l = excess_node(node -> left);
        Excess ex_r = excess_node(node -> right);
        
        Excess ex_l2 = excess_node_right(node -> left);
        Excess ex_r2 = excess_node_right(node -> right);
        ex_l.m_r = ex_l2.m_r;
        ex_r.m_r = ex_r2.m_r;
        node -> excess.combine(ex_l, ex_r);
    }

    public:
        DynamicTreeBP() {
            //start with only a leaf
            Leaf *l = new Leaf();
            root = TaggedPointer(l, IS_LEAF);
        }

        ~DynamicTreeBP() {
            auto [leaf_pt, node_pt, mark] = cast_pointer(root);
            if(mark == IS_LEAF) { delete leaf_pt; }
            else                { delete node_pt; }
        }
        
        // this methods is not used, since tree always starts from root
        DynamicTreeBP(std::vector<uint64_t> &bits, uint32_t num_bits) {
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
            recompute_excess(node);
            return {TaggedPointer(node, IS_NODE), l_num + r_num, l_ones + r_ones};
        }

        void only_root() {
            insert(0, false);
            insert(0, true);
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
                memory_total = memory_leafs + memory_nodes + sizeof(DynamicTreeBP);
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

        /* rank and select */
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

        /* rank and select */

        /* fwd_search*/
        int_type fwd_search(int_type i, int32_t d) {
            int_type bits_left = 0;
            /* traverse to block containing i */
            TaggedPointer cur = root;
            std::stack<std::pair<Node*, bool>> stack;
            bool left_right = false;
            while(is_node(cur)) {
                Node *node = get_node(cur);
                if(i - bits_left < node -> num) {
                    cur = node -> left;
                    left_right = LEFT;
                }
                else {
                    cur = node -> right;
                    left_right = RIGHT;
                    bits_left += node -> num;
                }
                stack.push({node, left_right});
            }
            /* traverse to block containing i */

            Leaf *leaf = get_leaf(cur);
            auto[excess, index_block] = fwd_search_block(leaf, i - bits_left, d);
            if(index_block != NOT_FOUND) {
                DBG << "found in first block " << V(excess) _ V(index_block) _ V(i + index_block) _ V(bits_left) _ V(d) _ V(i) << "\n";
                return bits_left + index_block;
            } 
            DBG << "did not found in first block " << V(excess) _ V(bits_left) _ V(d) _ V(i) << "\n";
            /* go up tree until found right subtree */
            while(!stack.empty()) {
                auto[node, lr] = stack.top(); stack.pop();
                if(lr == LEFT) {
                    Excess ex_right = excess_node(node -> right);
                    // std::cout << "left " << V(excess) _ V(ex_right.m) _ V(bits_left) << "\n";
                    assert(d != excess);
                    if((excess > d && excess + ex_right.m <= d) || (excess < d && excess + ex_right.m >= d)) {
                        cur = node -> right; //found subtree to descend
                        bits_left += node -> num;
                        break;
                    }
                    excess += ex_right.e;
                }
                else {
                    //came from right
                    bits_left -= node -> num;
                    // std::cout << "right " << V(excess) _ V(bits_left) << "\n";
                }
            }
            /* go up tree until found right subtree */
            
            /* descend into correct subtree subtree */
            while(is_node(cur)) {
                Node *node = get_node(cur);
                Excess ex_left = excess_node(node -> left);
                assert(d != excess); //otherwise would have found it already
                if((excess > d && excess + ex_left.m <= d) || (excess < d && excess + ex_left.m >= d)) { //go left
                    cur = node -> left;
                }
                else { //go right
                    excess += ex_left.e;
                    bits_left += node -> num;
                    cur = node -> right;
                }
            }
            /* descend into correct subtree subtree */
            Leaf *leaf2 = get_leaf(cur);
            // auto[d1, j] = fwd_search_block(leaf2, 0, -excess);
            auto[d1, j] = fwd_search_block(leaf2, 0, d - excess);
            DBG << "found in second block " << V(excess) _ V(bits_left) _ V(cur_size) _ V(j) _ V(d) _ V(i) << "\n";
            assert(j != NOT_FOUND);
            return bits_left + j;
        }


        std::pair<int32_t, int_type> fwd_search_block(Leaf *leaf, int_type i, int32_t d) {
            int32_t d1 = 0;
            for(int_type j = i; j < leaf -> size(); j++) {
                d1 += leaf -> access(j)? 1 : -1;
                if(d1 == d) return {d, j};
            }
            return {d1, NOT_FOUND};
        }
        /* fwd_search*/

        /* bwd_search*/
        //similar to fwd_search, but from right to left
        int_type bwd_search(int_type i, int32_t d) {
            int_type bits_left = 0;
            /* traverse to block containing i */
            TaggedPointer cur = root;
            std::stack<std::pair<Node*, bool>> stack;
            bool left_right = false;
            while(is_node(cur)) {
                Node *node = get_node(cur);
                if(i - bits_left < node -> num) {
                    cur = node -> left;
                    left_right = LEFT;
                }
                else {
                    cur = node -> right;
                    left_right = RIGHT;
                    bits_left += node -> num;
                }
                stack.push({node, left_right});
            }
            /* traverse to block containing i */

            Leaf *leaf = get_leaf(cur);
            auto[excess, index_block] = bwd_search_block(leaf, i - bits_left, d);
            if(index_block != NOT_FOUND) {
                // std::cout << "found in first block " << V(excess) _ V(index_block) _ V(i + index_block) _ V(bits_left) _ V(d) _ V(i) << "\n";
                return bits_left + index_block;
            } 
            // std::cout << "did not found in first block " << V(stack.size()) _ V(excess) _ V(bits_left) _ V(d) _ V(i) << "\n";
            
            /* go up tree until found right subtree */
            // std::cout << "going up tree \n";
            while(!stack.empty()) {
                auto[node, lr] = stack.top(); stack.pop();
                if(lr == RIGHT) {
                    // std::cout << "right up \n";
                    bits_left -= node -> num;
                    Excess ex_left = excess_node_right(node -> left);
                    assert(d != excess);
                    // std::cout << V(excess) _ V(ex_left.m_r) _ V(d) << "\n";
                    if((excess > d && excess + ex_left.m_r <= d) || (excess < d && excess + ex_left.m_r >= d)) {
                        // std::cout << "found left subtree \n";
                        cur = node -> left; //found subtree to descend
                        break;
                    }
                    excess += ex_left.e;
                }
                else {
                    //  std::cout << "left up \n";
                    //form left -> nothing to do
                }
            }
            /* go up tree until found right subtree */
            
            /* descend into correct subtree subtree */
            while(is_node(cur)) {
                Node *node = get_node(cur);
                Excess ex_right = excess_node_right(node -> right);
                assert(d != excess); //otherwise would have found it already
                // std::cout << V(excess) _ V(ex_right.m_r) _ V(d) << "\n";
                if((excess > d && excess + ex_right.m_r <= d) || (excess < d && excess + ex_right.m_r >= d)) { //go right
                    // std::cout << "descend right \n";
                    cur = node -> right;
                    bits_left += node -> num;
                }
                else { //go left
                    // std::cout << "descend left \n";
                    excess += ex_right.e;
                    cur = node -> left;
                }
            }
            /* descend into correct subtree subtree */
            Leaf *leaf2 = get_leaf(cur);
            // auto[d1, j] = bwd_search_block(leaf2, leaf2 -> size() - 1, -excess); //d - excess?
            auto[d1, j] = bwd_search_block(leaf2, leaf2 -> size() - 1, d - excess); //d - excess?
            // std::cout << "found in second block " << V(excess) _ V(bits_left) _ V(j) _ V(d) _ V(i) << "\n";
            assert(j != NOT_FOUND);
            return bits_left + j;
        }

        std::pair<int32_t, int_type> bwd_search_block(Leaf *leaf, int_type i, int32_t d) {
            int32_t d1 = 0;
            for(int32_t j = (int32_t)i; j >= 0; j--) {
                d1 -= leaf -> access(j)? 1 : -1; // -= !!!
                if(d1 == d) return {d, j};
            }
            return {d1, NOT_FOUND};
        }
        /* bwd_search*/



        /* tree operations*/

        int_type close(int_type v) {
            return fwd_search(v, 0);
        }

        int_type open(int_type v) {
            return bwd_search(v, 0);
        }

        int_type enclose(int_type v) {
            return bwd_search(v, -2);
        }

        int_type ith_node(int_type i) {
            return select(i, true);
        }

        int_type next_sibling(int_type v) {
            return fwd_search(v, 1);
        }

        int_type subtree_size(int_type v) {
            return (close(v) - v + 1) / 2;
        }

        int_type parent(int_type v) {
            assert(v != 0);
            return enclose(v);
        }

        int_type child_i(int_type v, int_type i) {
            int_type child = v + 1;
            int_type count = 0;
            while(access(child)) {
                count++;
                if(count == i) return child;
                child = close(child) + 1;
            }
            assert(false); //child does not exist
            return 0;
        }

        int_type children(int_type v) {
            assert(access(v));
            if(!access(v + 1)) return 0; //leaf
            int_type k = 1; //atleast one child
            int_type child = v + 1;
            int_type end_child = close(child);
            while(end_child + 1 < cur_size && access(end_child + 1)) {
                k++;
                child = end_child + 1;
                end_child = close(child);
            }
            return k;
        }

        //node v, pos i, pos j, num children
        std::tuple<int_type, int_type, int_type> child_i_j_num_child(int_type v, int_type i, int_type j) {
            assert(access(v)); //v is (
            int_type pos_i = 0, pos_j = 0, num_child = 0;
            if(!access(v + 1)) return {pos_i, pos_j, num_child}; //leaf
            num_child = 1;
            int_type child = v + 1;  //( of first child
            int_type end_child = close(child); // ) of first child
            if(i == 1) pos_i = child; 
            if(j == 1) pos_j = child; 
            //0 is never a child, since it is the root
            while((pos_i == 0 || pos_j == 0) && end_child + 1 < cur_size && access(end_child + 1)) {
                // std::cout << V(child) _ V(end_child) _ V(num_child) << "\n";
                num_child++;
                child = end_child + 1;
                end_child = close(child);
                if(i == num_child) pos_i = child; 
                if(j == num_child) pos_j = child; 
            }
            return {pos_i, pos_j, num_child};
        }

        void insert_child(int_type v, int_type i, int_type k) {
            // std::cout << "insert child " << V(v) _ V(i) _ V(k) << "\n";
            auto[l, r, num_childs] = child_i_j_num_child(v, i, i + k);
            // std::cout << V(l) _ V(r) _ V(num_childs) << "\n";
            l = (i <= num_childs)? l: close(v); // + 1?
            r = (i + k <= num_childs)? r : close(v);
            // std::cout << V(l) _ V(r) _ V(num_childs) << "\n";
            insert(r, 0); //right first
            insert(l, 1);
        }

        void delete_node(int_type v) {
            assert(access(v));
            assert(v != 0); //not root
            remove(close(v)); //right first
            remove(v);
        }
        /* tree operations*/

        //updates, height, balance and excess
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
                recompute_excess(node);
                stack.pop(); 
            }
            if(child != NULL && is_imbalanced(child)) {
                root = TaggedPointer(rebalance(child), IS_NODE); 
                if(is_node(root)) recompute_excess(get_node(root));
            }
        }

        void update_excess_path(std::stack<std::pair<Node*, bool>> &stack) {
            while(!stack.empty()) { 
                auto[node, path_dir] = stack.top(); stack.pop();
                recompute_excess(node); 
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
            if(leaf -> is_full()) {
                TaggedPointer pt = split_leaf(leaf);
                if(stack.size() == 0) { root = pt;} //root pointed to leaf, now to splitter
                else {
                    auto[parent, lr] = stack.top();
                    assign_pointer(pt, parent, lr);
                }
                update_path(stack);
            }
            else {
                update_excess_path(stack);
            }
        }

        TaggedPointer split_leaf(Leaf* leaf) {
            assert(leaf -> is_full());
            DBG << "split \n";
            Leaf *r_leaf = new Leaf();
            leaf -> split(*r_leaf);
            Node *splitter = new Node(leaf, r_leaf);
            recompute_excess(splitter);
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
                    recompute_excess(node);
                }
            };
            if(is_leaf(root)) {  //do nothing
                DBG << "only root \n";
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

            auto try_balance = [&](Leaf* a, Leaf* b) {
                StaticSmallBV<N>::simple_balance(*a, *b);
                return true;
            };

            auto merge_leaf = [&](Leaf* target, Leaf* other) {
                target -> simple_append(*other);
            };
            auto sum_small =  [&](Leaf* a, Leaf* b) { return a -> size() + b -> size() <= (MAX_SIZE / 2); };
            bool merged = false;
            //switched order: merge first, than try balance
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
                        if(grandparent != NULL) {
                            DBG << "grandparent gets leaf \n";
                            assign_pointer(pt_leaf, grandparent, lr_grandparent);
                            recompute_excess(grandparent);
                        }
                        else {
                            root = pt_leaf;
                            if(is_node(root)) {
                                recompute_excess(get_node(root));
                            }
                        }                   
                        parent -> invalidate_pointers();
                        delete parent;
                        delete leaf2;
                        merged = true;
                    }
                    else if(try_balance(leaf, leaf2)) {
                        DBG << "balance left a) \n";
                        update_stack.pop(); //delete parent from update_stack, num and ones got updated already
                        parent -> reinit_from_leaf(leaf);
                        recompute_excess(parent);
                    } 
                }
                else { // b)
                    Node* inner1 = get_node(pt);
                    Leaf* leaf2 = get_leaf(inner1 -> left);
                    if(sum_small(leaf, leaf2)) {
                        DBG << "merge left b) \n";
                        merge_leaf(leaf, leaf2);
                        parent -> right = TaggedPointer(get_leaf(inner1 -> right), IS_LEAF);
                        update_stack.pop(); //delete parent from update_stack, num and ones got updated already
                        parent -> reinit_from_leaf(leaf);
                        recompute_excess(parent);
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
                        recompute_excess(inner1);
                        recompute_excess(parent);
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
                        // leaf2 -> append_bv(*leaf); //merge
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
                        recompute_excess(parent);
                    }
                }
                else { // b)
                    Node* inner1 = get_node(pt);
                    Leaf* leaf2 = get_leaf(inner1 -> right);
                    Leaf* leaf3 = get_leaf(inner1 -> left);
                    if(sum_small(leaf2, leaf)) {
                        DBG << "merge right b) \n";
                        // leaf2 -> append_bv(*leaf); //merge
                        merge_leaf(leaf2, leaf);
                        parent -> left = TaggedPointer(leaf3, IS_LEAF);
                        parent -> right = TaggedPointer(leaf2, IS_LEAF);
                        parent -> reinit_from_leaf(leaf3);
                        recompute_excess(parent);
                        inner1 -> invalidate_pointers();
                        delete inner1;
                        delete leaf;
                        merged = true;
                    }
                    else if(try_balance(leaf2, leaf)) {
                        DBG << "balance right b) \n";
                        parent -> num = leaf3 -> size() + leaf2 -> size();
                        parent -> ones = leaf3 -> num_ones() + leaf2 -> num_ones();
                        recompute_excess(inner1);
                        recompute_excess(parent); //maybe redundant
                    } 
                }
            }
            update_num_ones();
            update_path(stack); //always update path to recompute excess
            if(merged) {/* */}
        }

        //returns new subtree root
        Node* rebalance(Node* z) { 
            DBG << "rebalancing \n";
            // std::cout << "rebalancing \n";
            // pr_tree();
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

        uint64_t size() {
            return cur_size;
        }

        void pr_tree() {
            pr_tree_rec(root, 0);
        }

        //|balance| <= 1
        bool check_invariant() {
            return check_invariant_rec(root);
        }

        bool check_node_excess() {
            return check_node_excess_rec(root);
        }

        uint32_t blocks_per_leaf() {return N;}

        std::vector<uint32_t> all_out_deg() {
            std::vector<uint32_t> out_deg(cur_size / 2, 0);
            std::stack<uint32_t> node_stack;
            uint32_t counter = 0;
            rec_all_out_deg(root, out_deg, node_stack, counter);
            return out_deg;
        }

        void rec_all_out_deg(TaggedPointer pt, std::vector<uint32_t> &out_deg, std::stack<uint32_t> &node_stack, uint32_t &counter) {
            auto[leaf, node, mark] = cast_pointer(pt);
            if(mark == IS_LEAF) {
                //dfs traversal, since bp is in dfs order
                for(uint32_t j = 0; j < leaf -> size(); j++) {
                    if(leaf -> access(j)) { // ( -> traverse down in tree
                        node_stack.push(counter); //node id
                        counter++;
                    }
                    else { // ) -> traverse up in tree
                        node_stack.pop();
                        if(!node_stack.empty()) { //not root
                            out_deg[node_stack.top()]++; //increase out degree of parent
                        }
                    }
                }
            }
            else {
                rec_all_out_deg(node -> left, out_deg, node_stack, counter);
                rec_all_out_deg(node -> right, out_deg, node_stack, counter);
            }
        }

    uint64_t leaf_blocks() { return N;}

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

        bool check_node_excess_rec(TaggedPointer pt) {
            auto[leaf, node, mark] = cast_pointer(pt);
            if(mark == IS_LEAF) return true;
            if(!check_node_excess_rec(node -> left)) return false;
            if(!check_node_excess_rec(node -> right)) return false;

            Excess node_ex = node -> excess;

            Excess ex_l = excess_node(node -> left);
            Excess ex_r = excess_node(node -> right);
            Excess ex_l2 = excess_node_right(node -> left);
            Excess ex_r2 = excess_node_right(node -> right);
            ex_l.m_r = ex_l2.m_r;
            ex_r.m_r = ex_r2.m_r;
            Excess comb;
            comb.combine(ex_l, ex_r);
            bool ok = true;
            ok &= comb.e == node_ex.e;
            if(!ok) std::cout << V(comb.e) _ V(node_ex.e) << "\n";
            ok &= comb.m == node_ex.m;
            if(!ok) std::cout << V(comb.m) _ V(node_ex.m) << "\n";
            ok &= comb.m_r == node_ex.m_r;
            if(!ok) std::cout << V(comb.m_r) _ V(node_ex.m_r) << "\n";
            return ok;
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

            recompute_excess(y);
            recompute_excess(x);
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

            recompute_excess(x);
            recompute_excess(y);
            return y; //new root
        }

        void pr_tree_rec(TaggedPointer pt, int_type level) {
            auto [leaf_pt, node_pt, mark] = cast_pointer(pt);
            if(mark == IS_LEAF) {
                // std::cout << std::string(level + 1, ' ' ) << "[" << leaf_pt -> size() << "]\n";
                Excess ex_leaf = excess_leaf(*leaf_pt);
                Excess ex_leaf_right = excess_leaf_right(*leaf_pt);
                std::cout << std::string(level + 1, ' ' ) << "[" << leaf_pt -> size() << "] " 
                << "e=" << ex_leaf.e << " m=" << ex_leaf.m  << " m_r=" << ex_leaf_right.m_r
                << "\n";
                leaf_pt -> pr_all();
                // excess_leaf(*leaf_pt).print();
                return;
            }
            std::cout << std::string(level, ' ' ) << "( height, num, ones:" 
            << " " << node_pt ->height 
            << " " << node_pt ->num 
            << " " << node_pt ->ones 
            << " e=" << node_pt -> excess.e
            << " m=" << node_pt -> excess.m
            << " m_r=" << node_pt -> excess.m_r
            << "\n"; 
            // node_pt -> excess.print();
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

        const static int_type NOT_FOUND = std::numeric_limits<int_type>::max();;
};
