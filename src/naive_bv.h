#pragma once

#include <vector>
#include <cassert>
#include <cstdint>
#include <bitset>

#include "src/utility.h"

/* Naive bitvector implementation for checking correctness. */
class NaiveBV {
    public:
        using int_type = uint32_t;

        NaiveBV() {};

        NaiveBV(std::vector<bool> &v) : arr(v) {};

        bool access(int_type i) {
            assert(i < arr.size());
            return arr[i];
        }

        void insert(int_type i, bool bit) {
            assert(i < arr.size() + 1);
            arr.insert(arr.begin() + i, bit);
            return;
        }

        void remove(int_type i) {
            assert(i < arr.size());
            arr.erase(arr.begin() + i);
            return;
        }

        void flip(int_type i) {
            assert(i < arr.size());
            arr[i] = arr[i] ^ 1;
            return;
        }

        // b == true -> rank 1, b == false -> rank 0
        int_type rank(int_type i, bool b) {
            assert(i < arr.size());
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

        int_type size() {
            return arr.size();
        }

        void print() {
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
