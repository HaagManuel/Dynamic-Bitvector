#pragma once

#include <cstdint>

#include "src/oop/node.h"

template<unsigned int N = 2> 
class Leaf final : public Node {

    public:
        using int_type = uint32_t;

        Leaf() : bv() {};

        Leaf(StaticSmallBV<N> &&bv_in) : bv(bv_in) {}

        ~Leaf() {
            bv.~StaticSmallBV<N>();
        }
        bool access(int_type i) override { return bv.access(i); }
        bool flip(int_type i) override { bv.flip(i); return bv.access(i);}
        Node* insert(int_type i, bool b) override { 
            bv.insert(i, b); 
            if(bv.is_full()) {
                StaticSmallBV<N> new_bv;
                bv.split(new_bv);
                Leaf* r_leaf = new Leaf<N>(std::move(new_bv));
                return new InnerNode(this, r_leaf);
            }
            return this;
        }
        
        int_type rank1(int_type i, int_type ones_before) override { return ones_before + bv.rank(i, true); }
        int_type select1(int_type i, int_type bits_before) { return bits_before + bv.select(i, true); }
        int_type select0(int_type i, int_type bits_before) { return bits_before + bv.select(i, false); }

        int_type get_num() override { return bv.size();} 
        int_type get_ones() override { return bv.rank(bv.size(), true);} 
        int_type get_height() override { return 0;} 
        int32_t get_balance() override { return 0;} 
        void update_height() override { /* do nothing */} 
        void set_ones(int_type) override { /* do nothing */} 
        void set_num(int_type) override { /* do nothing */} 

        void pr_tree(int_type level) override { 
            std::cout << std::string(level + 1, ' ' ) << "[" << bv.size() << "]\n";
        } 
        bool check_invariant() override {return true;}

    private:
        StaticSmallBV<N> bv;
};