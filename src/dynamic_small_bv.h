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

#define DBG if constexpr(DEBUG_OUT) std::cout
#define ONLY_DBG if constexpr(DEBUG_OUT)

/* Small BV that reallocates array instead of using a static size.
   More space efficient, but slower. Implementaiton not finally tested.
*/
//N = max number of blocks -> than other algorithm can call split
template<unsigned int N = 4> 
class DynamicSmallBV {
    private:
        static const bool DEBUG_OUT = false; 
    public:
        DynamicSmallBV() {
            arr = std::vector<uint64_t>(MAX_FREE, 0);
            cur_size = 0;
        }

        DynamicSmallBV(uint32_t n) {
            arr = std::vector<uint64_t>(n, 0);
            cur_size = 0;
        }

        void increase_memory(uint32_t n) {
            //dont allocate memory, if max size is reached, leaf is split before in tree
            n = std::min(n, (uint32_t) (N - arr.capacity()));
            if(n == 0) return;
            DBG << "increasing memory by " << V(n);
            DBG << V(arr.capacity()) << " ";
            arr.reserve(arr.size() + n);
            arr.resize(arr.size() + n);
            DBG << V(arr.capacity()) << "\n";
        }

        void decrease_memory(uint32_t n) {
            assert(n > 0);
            DBG << "decreasing memory by " << V(n);
            DBG << V(arr.capacity()) << " ";
            //reallocate smaller vector an copy contents
            std::vector<uint64_t>(arr.begin(), std::max(arr.begin() + 1, arr.end() - n)).swap(arr);
            DBG << V(arr.capacity()) << "\n";
        }

        void set_memory_to(uint32_t n) {
            if(n > arr.size()) {
                increase_memory(n - arr.size());
            }
            else if(n < arr.size()) {
                decrease_memory(arr.size() - n);
            }
        }

        //access each block from right to left
        bool access(uint32_t i) {
            auto[block, pos] = get_block_and_index(i);
            // std::cout << V(i) _ V(block) _ V(arr.size()) _ V(pos) _ V(BITS) << "\n";
            assert(block < arr.size() && pos < BITS);
            return (arr[block] >> pos) & (1ULL);
        }

        void insert(uint32_t i, bool bit) {
            assert(i < cur_size + 1); //allowed to insert direclty behind last
            auto[block, pos] = get_block_and_index(i);
            uint32_t last_touched_block = cur_size / BITS;
            cur_size++;
            if(get_cur_block() >= arr.size()) increase_memory(MAX_FREE);

            assert(block < arr.size());

            //save last bit of last touched block, except it is the last block -> there is always a next block for dynamic bv
            if(last_touched_block != block) {
                // if(last_touched_block + 1 != cur_blocks) { //not last block
                //if last touched block is full, need to safe last bit
                    arr[last_touched_block + 1] |= (arr[last_touched_block] >> (BITS - 1)); 
                // }
                arr[last_touched_block] <<= 1; //left shift, since we access from right to left
            }

            //shift cur_blocks right of target block
            //int32_t -> avoid underflow
            for(int32_t j = (int32_t)last_touched_block - 1; j >= (int32_t)block + 1; j--) {
                arr[j + 1] |= (arr[j] >> (BITS - 1)); //save bit from before
                arr[j] <<= 1;
            }
            // if(block + 1 != cur_blocks) {
            //     arr[block + 1] |= (arr[block] >> (BITS - 1));
            // }
            arr[block + 1] |= (arr[block] >> (BITS - 1));
            
            //shift part of target block 
            uint64_t right_part = arr[block] & ((1ULL << pos) - 1);
            arr[block] <<= 1;
            arr[block] &= ~(1ULL << pos); //clear place for new bit
            arr[block] |= ((uint64_t)bit << pos); //inserted bit
            arr[block] &= ~((1ULL << pos) - 1); //clear right part
            arr[block] |= right_part; //restore right part
            
            return;
        }

        void remove(uint32_t i) {
            assert(i < cur_size && cur_size > 0);
            auto[block, pos] = get_block_and_index(i);

            uint32_t last_touched_block = cur_size / BITS; 
            // uint32_t cur_blocks = num_blocks();
            if(cur_size == MAX_SIZE) last_touched_block--; //when full

            assert(last_touched_block < num_blocks());

            //target block
            uint64_t right_part = arr[block] & ((1ULL << pos) - 1);
            arr[block] >>= 1; //read bits from right to left
            arr[block] &= ~((1ULL << pos) - 1); //clear right part
            arr[block] |= right_part; //restore right part

            //shift cur_blocks right of target block
            //save first bit from right block to last bit (0 after shift) of current block
            for(uint32_t j = block; j < last_touched_block; j++) {
                // assert(j < cur_blocks);
                arr[j] |= ((arr[j + 1] & 1ULL) << (BITS - 1)); //take first bit, set last bit
                arr[j + 1] >>= 1;
            }
            cur_size--;
            if(num_free_blocks() > MAX_FREE) {
                decrease_memory(MAX_FREE);
            }
            return;
        }

        void flip(uint32_t i) {
            auto[block, pos] = get_block_and_index(i);
            assert(block < arr.size());
            arr[block] ^= (1ULL << pos);
            return;
        }

        uint32_t rank(uint32_t i, bool b) { 
            auto[block, pos] = get_block_and_index(i);
            int32_t count1 = 0;
            for(uint32_t j = 0; j < block; j++) {
                count1 += std::popcount(arr[j]);
            }
            //rank excluding pos (shift 64 is undefined)
            if(pos > 0) count1 += std::popcount((arr[block] << (BITS - pos))); //without bits right of pos
            return b * count1 + (b^1) * (i - count1); //rank1 or rank0
        }

        uint32_t select(uint32_t i, bool b) {
            assert(i > 0); //beginn at 1.
            return b? select1(i) : select0(i);
        }

        uint32_t size() {
            return cur_size;
        }

        uint32_t max_size() {
            return MAX_SIZE;
        }

        bool has_capacity() {
            return cur_size < MAX_SIZE;
        }

        bool is_full() { return cur_size == MAX_SIZE; }

        bool is_empty() { return cur_size == 0; }

        bool is_aligned() { return (cur_size & (BITS - 1)) == 0;}

        uint64_t get_word(uint32_t block_id) {
            assert(block_id < arr.size());
            return arr[block_id];
        }

        //size must be set accordingly after
        void set_word(uint64_t w, uint32_t block_id) {
            assert(block_id < arr.size());
            arr[block_id] = w;
        } 

        void or_word(uint64_t w, uint32_t block_id) {
            assert(block_id < arr.size());
            arr[block_id] |= w;
        }

        void shift_block_left(uint32_t block_id, uint32_t k) {
            assert(k < BITS && block_id < arr.size());
            arr[block_id] <<= k;
        }

        void shift_block_right(uint32_t block_id, uint32_t k) {
            assert(k < BITS && block_id < arr.size());
            arr[block_id] >>= k;
        }

        

        //bits in new range are not cleared!
        void set_size(uint32_t new_size) {
            assert(new_size <= MAX_SIZE);
            cur_size = new_size;
            set_memory_to(get_cur_block() + 1);
        }

        void add_size(int32_t add) {
            assert(cur_size + add < MAX_SIZE);
            cur_size += add;
        }
        //copy right half of this vector to other, other should have the right size
        void split(DynamicSmallBV<N> &other) {
            assert(cur_size == MAX_SIZE); //otherwise not intended to split
            assert(N % 2 == 0); //for even split
            uint32_t j = 0;
            other.set_size(cur_size / 2);
            for(uint32_t block = N / 2; block < N; block++) {
                other.set_word(arr[block], j);
                j++;
            }
            set_size(cur_size / 2); //throw away half of bits, not cleared
            // cur_size /= 2; 
            // other.set_size(cur_size);
        }

        uint32_t num_ones() {
            return rank(cur_size, true);
        }


        //simpler version for correctness but slow

        //appends other to this in a naive way
        void simple_append(DynamicSmallBV<N> &other) {
            assert(cur_size + other.size() <= MAX_SIZE);
            for(uint32_t i = 0; i < other.size(); i++) {
                insert(size(), other.access(i));
            }
        }
        //balanaced bv's in naive way
        static void simple_balance(DynamicSmallBV<N> &left, DynamicSmallBV<N> &right) {
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

        void pr_all() {
            for(uint32_t i = 0; i <= last_used_block(); i++) {
                utility::pr_bit(arr[i]);
            }
        }


        uint32_t last_used_block() {
            uint32_t cur_blocks = cur_size / BITS;
            if(is_aligned()) cur_blocks--;
            return cur_blocks;
        }
        uint32_t get_cur_block() {return cur_size / BITS;}
        uint32_t num_blocks() { return arr.size(); }
        uint32_t num_free_blocks() { return num_blocks() - get_cur_block(); }


        //dummy method
        static void balance_right_to_left(DynamicSmallBV<N>, DynamicSmallBV<N>) { return; }

        static void balance_left_to_right(DynamicSmallBV<N>, DynamicSmallBV<N>) { return; }

        void append_bv(DynamicSmallBV<N>) { return; }

    static const uint32_t BITS = 64;
    static const uint32_t MAX_SIZE = N * BITS; 

    private:

        inline const std::pair<uint32_t, uint32_t> get_block_and_index(uint32_t i) {
            return {i / BITS, i & (BITS - 1)};
        }

        uint32_t select1(uint32_t i) {
            uint32_t count1 = 0;
            for(uint32_t j = 0; j <= get_cur_block(); j++) {
                uint32_t next = std::popcount(arr[j]);
                if(count1 + next >= i) { //found right block
                    for(uint32_t k = 0; k < BITS; k++) { //scan for index
                        uint32_t index = j * BITS + k;
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

        uint32_t select0(uint32_t i) {
            uint32_t count0 = 0;
            for(uint32_t j = 0; j <= get_cur_block(); j++) {
                uint32_t next = std::popcount(~arr[j]);
                if(count0 + next >= i) { //found right block
                    for(uint32_t k = 0; k < BITS; k++) { //scan for index
                        uint32_t index = j * BITS + k;
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
        
        std::vector<uint64_t> arr; //size() of vector -> num blocks
        uint32_t cur_size;
        
        static const uint32_t MAX_FREE = 1;
        // static const uint32_t MAX_FREE = 4;
};
