#pragma once

#include <vector>
#include <cassert>
#include <cstdint>
#include <array>
#include <utility>
#include <iostream>
#include <bitset>
#include <bit>

#include "src/utility.h"

/* 
 Small bitvector of fixed size.
 Implements leaf operations of dynamic bitvector.
 Holds N uint64_t.
*/

template<unsigned int N = 4> 
class StaticSmallBV {
    public:
        using int_type = uint32_t;

        StaticSmallBV() {
            std::fill(arr.begin(), arr.end(), 0);
        }

        //access each block from right to left
        bool access(int_type i) {
            assert(i < MAX_SIZE);
            auto[block, pos] = get_block_and_index(i);
            return (arr[block] >> pos) & (1ULL);
        }

        void insert(int_type i, bool bit) {
            assert(i < cur_size + 1); //allowed to insert direclty behind last
            assert(cur_size < MAX_SIZE);
            auto[block, pos] = get_block_and_index(i);

            int_type last_touched_block = cur_size / BITS;
            assert(last_touched_block < BLOCKS);

            //save last bit of last touched block, except it is the last block
            if(last_touched_block != block) {
                if(last_touched_block + 1 != BLOCKS) { //not last block
                //if last touched block is full, need to safe last bit
                    arr[last_touched_block + 1] |= (arr[last_touched_block] >> (BITS - 1)); 
                }
                arr[last_touched_block] <<= 1; //left shift, since we access from right to left
            }

            //shift blocks right of target block
            //int32_t -> avoid underflow
            for(int32_t j = (int32_t)last_touched_block - 1; j >= (int32_t)block + 1; j--) {
                assert(j >= 0 && j < (int32_t) BLOCKS);
                arr[j + 1] |= (arr[j] >> (BITS - 1)); //save bit from before
                arr[j] <<= 1;
            }
            if(block + 1 != BLOCKS) {
                arr[block + 1] |= (arr[block] >> (BITS - 1));
            }
            //shift part of target block 
            uint64_t right_part = arr[block] & ((1ULL << pos) - 1);
            arr[block] <<= 1;
            arr[block] &= ~(1ULL << pos); //clear place for new bit
            arr[block] |= ((uint64_t)bit << pos); //inserted bit
            arr[block] &= ~((1ULL << pos) - 1); //clear right part
            arr[block] |= right_part; //restore right part
            cur_size++;
            return;
        }

        void remove(int_type i) {
            assert(i < cur_size && cur_size > 0);
            auto[block, pos] = get_block_and_index(i);

            int_type last_touched_block = cur_size / BITS; 
            if(cur_size == MAX_SIZE) last_touched_block--; //when full
            assert(last_touched_block < BLOCKS);

            //target block
            uint64_t right_part = arr[block] & ((1ULL << pos) - 1);
            arr[block] >>= 1; //read bits from right to left
            arr[block] &= ~((1ULL << pos) - 1); //clear right part
            arr[block] |= right_part; //restore right part

            //shift blocks right of target block
            //save first bit from right block to last bit (0 after shift) of current block
            for(int_type j = block; j < last_touched_block; j++) {
                assert(j < BLOCKS);
                arr[j] |= ((arr[j + 1] & 1ULL) << (BITS - 1)); //take first bit, set last bit
                arr[j + 1] >>= 1;
            }
            cur_size--;
            return;
        }

        void flip(int_type i) {
            auto[block, pos] = get_block_and_index(i);
            arr[block] ^= (1ULL << pos);
            return;
        }

        // b == true -> rank 1, b == false -> rank 0
        int_type rank(int_type i, bool b) {
            if(i == MAX_SIZE) {
                int32_t count1 = 0;
                for(int_type j = 0; j < BLOCKS; j++) count1 += std::popcount(arr[j]);
                return count1;
            }
            auto[block, pos] = get_block_and_index(i);
            int32_t count1 = 0;
            for(int_type j = 0; j < block; j++) {
                count1 += std::popcount(arr[j]);
            }
            //rank excluding pos (shift 64 is undefined)
            if(pos > 0) count1 += std::popcount((arr[block] << (BITS - pos))); //without bits right of pos
            return b * count1 + (b^1) * (i - count1); //rank1 or rank0
        }

        int_type select(int_type i, bool b) {
            assert(i > 0); //beginn at 1.
            return b? select1(i) : select0(i);
        }

        int_type size() {
            return cur_size;
        }

        int_type max_size() {
            return MAX_SIZE;
        }

        bool has_capacity() {
            return cur_size < MAX_SIZE;
        }

        bool is_full() { return cur_size == MAX_SIZE; }

        bool is_empty() { return cur_size == 0; }

        bool is_aligned() { return (cur_size & (BITS - 1)) == 0;}

        uint64_t get_word(int_type block_id) {
            assert(block_id < BLOCKS);
            return arr[block_id];
        }

        //size must be set accordingly after
        void set_word(uint64_t w, int_type block_id) {
            assert(block_id < BLOCKS);
            arr[block_id] = w;
        } 

        void or_word(uint64_t w, int_type block_id) {
            assert(block_id < BLOCKS);
            arr[block_id] |= w;
        }

        void shift_block_left(int_type block_id, int_type k) {
            assert(k < BITS && block_id < BLOCKS);
            arr[block_id] <<= k;
        }

        void shift_block_right(int_type block_id, int_type k) {
            assert(k < BITS && block_id < BLOCKS);
            arr[block_id] >>= k;
        }

        int_type last_used_block() {
            int_type blocks = cur_size / BITS;
            if(is_aligned()) blocks--;
            return blocks;
        }

        //bits in new range are not cleared!
        void set_size(int_type new_size) {
            assert(new_size <= MAX_SIZE);
            cur_size = new_size;
        }

        void add_size(int32_t add) {
            assert(cur_size + add < MAX_SIZE);
            cur_size += add;
        }
        //copy right half of this vector to other
        void split(StaticSmallBV<N> &other) {
            assert(cur_size == MAX_SIZE); //otherwise not intended to split
            assert(N % 2 == 0); //for even split
            int_type j = 0;
            for(int_type block = BLOCKS / 2; block < BLOCKS; block++) {
                other.set_word(arr[block], j);
                j++;
            }
            cur_size /= 2; //throw away half of bits, not cleared
            other.set_size(cur_size);
        }

        int_type num_ones() {
            return rank(cur_size, true);
        }


        //simpler version for correctness but slow

        //appends other to this in a naive way
        void simple_append(StaticSmallBV<N> &other) {
            assert(cur_size + other.size() <= MAX_SIZE);
            for(int_type i = 0; i < other.size(); i++) {
                insert(size(), other.access(i));
            }
        }
        //balanaced bv's in naive way
        static void simple_balance(StaticSmallBV<N> &left, StaticSmallBV<N> &right) {
            bool first = left.size() + 2 < right.size();
            while(left.size() + 2 < right.size()) {
                left.insert(left.size(), right.access(0));
                right.remove(0);
            }
            while(!first && left.size() > right.size() + 26) {
                right.insert(0, left.access(left.size() - 1));
                left.remove(left.size() - 1);
            }
        }


        //append bits of other bv to this bv, this bv must have a full subblock for efficient word copy
        void append_bv(StaticSmallBV<N> &other) {
            assert(cur_size + other.size() <= MAX_SIZE);
            align_left(*this, other);
            assert(is_aligned()); //to make use of word copy
            int_type free_block = cur_size / BITS;
            int_type j = 0;
            for(int_type i = 0; i < other.size(); i += BITS) {
                set_word(other.get_word(j), free_block);
                j++;
                free_block++;
            }
            cur_size += other.size();
        }
        
        //copy trailing bits to left than call aligned routing for balancing
        //does not work completly
        static void align_left(StaticSmallBV<N> &left, StaticSmallBV<N> &right) {
            if(left.is_aligned()) return;
            int_type last_left = left.last_used_block();
            int_type last_right = right.last_used_block();
            int32_t k = left.size() & (left.BITS - 1);
            assert(right.size() + k <= right.max_size());
            if(last_right != right.BLOCKS - 1)  {//not last block -> save last k bits
                uint64_t k_last = right.get_word(last_right) >> (left.BITS - k);
                right.set_word(k_last, last_right + 1);
            }
            right.shift_block_left(last_right, k);
            for(int32_t j = (int32_t)last_right - 1; j >= 0; j--) {
                uint64_t k_last = right.get_word(j) >> (left.BITS - k);
                right.or_word(k_last, j + 1);
                right.shift_block_left(j, k);
            }
            uint64_t k_last_left = left.get_word(last_left); //test, mask only first k bits
            right.or_word(k_last_left , 0);
            left.add_size(-k);
            right.add_size(k);
        }

        //right is small, left is aligned than right is copied left
        static void balance_right_to_left(StaticSmallBV<N> &left, StaticSmallBV<N> &right) {
            // align_left(left, right);
            assert(left.is_aligned()); //to make use of word copy
            int_type left_blocks = left.last_used_block() + 1, right_blocks = right.last_used_block() + 1;
            int_type total_blocks = left_blocks + right_blocks;
            if(left.size() < right.size() && right_blocks > 1) {
                int_type free_block = left.size() / left.BITS;
                int_type to_copy = (total_blocks / 2) - left_blocks;
                for(int_type i = 0; i < to_copy; i++) {
                    left.set_word(right.get_word(i), free_block + i); //transfer blocks from right to left
                }
                for(int_type i = to_copy; i < right_blocks; i++) {
                    right.set_word(right.get_word(i), i - to_copy); //shift blocks to the left
                }
                left.add_size(left.BITS * to_copy);
                right.add_size( (-1) * (int32_t)(left.BITS * to_copy));
            }
        }

        //left is small and  must be aligned
        static void balance_left_to_right(StaticSmallBV<N> &left, StaticSmallBV<N> &right) {
            align_left(left, right);

            int_type left_blocks = left.last_used_block() + 1, right_blocks = right.last_used_block() + 1;
            int_type total_blocks = left_blocks + right_blocks;
            if(left.size() > right.size() && left_blocks > 1) {
                assert(left.is_aligned()); //to make use of word copy
                int_type to_copy = (total_blocks / 2) - right_blocks;
                // std::cout << V(to_copy) _ V(left_blocks) _ V(right_blocks) << "\n";
                for(int32_t i = (int32_t)right_blocks - 1; i >= 0; i--) { //shift right blocks to the right
                    right.set_word(right.get_word(i), to_copy + i); 
                }
                for(int_type i = 0; i < to_copy; i++) { //transfer blocks from left to right
                    right.set_word(left.get_word(left_blocks + i - to_copy), i); 
                }
                left.set_size(left.size() - left.BITS * to_copy);
                right.set_size(right.size() + left.BITS * to_copy);
            }
            // std::cout << V(left.size()) _ V(right.size()) << "\n";
            // left.pr_all();
            // right.pr_all();
        }

        void pr_all() {
            for(uint32_t i = 0; i <= last_used_block(); i++) {
                utility::pr_bit(arr[i]);
            }
        }

    static const int_type BITS = 64;
    static const int_type BLOCKS = N;
    static const int_type MAX_SIZE = BITS * BLOCKS;

    private:
        inline const std::pair<int_type, int_type> get_block_and_index(int_type i) {
            assert(i < MAX_SIZE);
            assert((i / BITS) < BLOCKS);
            assert((i & (BITS - 1)) < BITS);
            return {i / BITS, i & (BITS - 1)};
        }

        int_type select1(int_type i) {
            int_type count1 = 0;
            for(int_type j = 0; j < BLOCKS; j++) {
                int_type next = std::popcount(arr[j]);
                if(count1 + next >= i) { //found right block
                    for(int_type k = 0; k < BITS; k++) { //scan for index
                        int_type index = j * BITS + k;
                        if(index >= cur_size) return MAX_SIZE; //no inserted bits anymore
                        if((arr[j] >> k) & 1ULL) {
                            count1++;
                            if(count1 == i) return index;
                        }
                    }
                }
                count1 += next;
            }
            return MAX_SIZE; //not found
        }

        int_type select0(int_type i) {
            int_type count0 = 0;
            for(int_type j = 0; j < BLOCKS; j++) {
                int_type next = std::popcount(~arr[j]);
                if(count0 + next >= i) { //found right block
                    for(int_type k = 0; k < BITS; k++) { //scan for index
                        int_type index = j * BITS + k;
                        if(index >= cur_size) return MAX_SIZE; //no inserted bits anymore
                        if(~(arr[j] >> k) & 1ULL) {
                            count0++;
                            if(count0 == i) return index;
                        }
                    }
                }
                count0 += next;
            }
            return MAX_SIZE; //not found
        }
        
        std::array<uint64_t, N> arr;
        int_type cur_size = 0;
};
