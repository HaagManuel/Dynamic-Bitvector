#pragma once

#include <vector>
#include <cassert>
#include <cstdint>
#include <bitset>

#include "src/utility.h"
/* Naive bp-tree implementation for checking correctness */
class NaiveTreeBP {
    public:
        using int_type = int32_t;

        NaiveTreeBP() {};

        NaiveTreeBP(std::vector<bool> &v) : arr(v) {};


        void only_root() {
            arr.resize(2);
            arr[0] = true;
            arr[1] = false;
        }

        //v must be ( to node
        //i is index in bitarray

        bool access(int_type i) {
            assert(i < size());
            return arr[i];
        }

        void insert(int_type i, bool bit) {
            assert(i < size() + 1);
            arr.insert(arr.begin() + i, bit);
            return;
        }

        void remove(int_type i) {
            assert(i < size());
            arr.erase(arr.begin() + i);
            return;
        }

        int_type rank(int_type i, bool b) {
            assert(i < size());
            return std::count(arr.begin(), arr.begin() + i, b);
        }

        int_type select(int_type i, bool b) {
            int_type count = 0, index = 0;
            for(bool x : arr) {
                if(x == b) count++;
                if(count == i) return index;
                index++;
            }
            return arr.size(); //not found
        }

        int_type excess(int_type i) {
            int_type excess = 0; //rank 1 - rank 0 = #( - #)
            for(int_type j = 0; j <= i; j++) {
                assert(excess >= 0); //parentheses expression
                excess += access(j) ? 1 : -1;
            }
            return excess;
        }

        int_type fwd_search(int_type i, int_type d) {
            assert(i < size());
            int_type d1 = 0;
            for(int_type j = i; j < size(); j++) {
                d1 += access(j) ? 1 : -1;
                if(d1 == d) return j;
            }
            return size();
        }

        int_type bwd_search(int_type i, int_type d) {
            assert(i < size());
            int_type d1 = 0;
            for(int_type j = i; j >= 0; j--) {
                d1 -= access(j) ? 1 : -1;
                if(d1 == d) return j;
            }
            return size();
        }

        int_type close(int_type v) {
            return fwd_search(v, 0);
        }

        int_type open(int_type v) {
            return bwd_search(v, 0);
        }
        
        int_type enclose(int_type v) {
            return bwd_search(v, -2);
        }

        //min excess(i ... j)
        int_type min_excess(int_type i, int_type j) {
            int_type ex = excess(i);
            int_type min_ex = ex;
            for(int32_t k = i + 1; k <= j; k++) {
                ex += access(k) ? 1 : -1;
                min_ex = std::min(min_ex, ex);
            }
            return min_ex;
        }

        //position of the t-th minimum in excess(i ... j)
        int_type min_select(int_type i, int_type j, int_type t) {
            int_type min_ex = min_excess(i, j);
            int_type ex = excess(i - 1);
            int_type count = 0; //first child starts at 1
            for(int32_t k = i; k <= j; k++) {
                ex += access(k) ? 1 : -1;
                if(ex == min_ex) count++;
                std::cout << V(k) _ V(ex) _ V(count) _ V(min_ex) << "\n";
                if(count == t) return k;
            }
            return j;
        }

        //num minimum in excess(i ... j)
        int_type min_count(int_type i, int_type j) {
            int_type min_ex = min_excess(i, j);
            int_type ex = excess(i - 1);
            int_type num = 0;
            for(int32_t k = i; k <= j; k++) {
                ex += access(k) ? 1 : -1;
                if(ex == min_ex) num++;
            }
            return num;
        }

        // int_type child_i(int_type v, int_type t) {
        //     return min_select(v, close(v) - 0, t);
        // }

        int_type child_i(int_type v, int_type i) {
            int_type child = v + 1;
            int_type count = 0;
            while(access(child)) {
                count++;
                if(count == i) return child;
                child = close(child) + 1;
            }
            return 0;
        }

        int_type parent(int_type v) {
            assert(v != 0);
            return enclose(v);
        }

        int_type children(int_type v) {
            return min_count(v, close(v) - 2);
        }

        void insert_child(int_type v, int_type i, int_type k) {
            int_type l, r;
            l = (i <= children(v))? child_i(v, i) + 0: close(v);
            r = (i + k <= children(v))? child_i(v, i + k) + 0  : close(v);
            // std::cout << V(l) _ V(r) _ V(children(v)) << "\n";
            insert(r, 0); //right first
            insert(l, 1);
        }

        void delete_node(int_type v) {
            assert(access(v));
            remove(close(v)); //right first
            remove(v);
        }

        bool is_parentheses_expr() {
            int_type ex = 0;
            bool ok = true;
            for(int_type j = 0; j < size(); j++) {
                ex += access(j) ? 1 : -1;
                ok &= (ex >= 0);
            }
            return ok;
        }

        int_type subtree_size(int_type v) {
            return (close(v) - v + 1) / 2;
        }

        int_type size() {
            return arr.size();
        }

        void print() {
            for(bool b : arr) {
                std::cout << b;
            }
            std::cout << "\n";
        }

        void print_bit() {
            uint64_t i = 0;
            uint64_t x = 0;
            for(bool b : arr) {
                x |= (1LL << i) * b;
                i++;
                if(i % 64 == 0) {
                    utility::pr_bit(x);
                    x = 0;
                    i = 0;
                }
            }
            if(i > 0) utility::pr_bit(x);
        }

    private:
        std::vector<bool> arr;
};
