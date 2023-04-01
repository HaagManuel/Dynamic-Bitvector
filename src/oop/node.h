#pragma once

#include <cstdint>

#include "src/static_small_bv.h"

class Node {
    public:
        using int_type = uint32_t;
        
        virtual ~Node(){} //make the base class destructor virtual

        virtual bool access(int_type i) = 0;
        virtual bool flip(int_type i) = 0; //returns if bit was set to true or false
        virtual Node* insert(int_type i, bool b) = 0;
        //remove
        virtual int_type rank1(int_type i, int_type ones_before) = 0;
        virtual int_type select1(int_type i, int_type bits_before) = 0;
        virtual int_type select0(int_type i, int_type bits_before) = 0;

        virtual int_type get_num() = 0;
        virtual int_type get_ones() = 0;
        virtual int_type get_height() = 0;
        virtual int32_t get_balance() = 0;
        virtual void update_height() = 0;
        virtual void set_num(int_type n) = 0;
        virtual void set_ones(int_type o) = 0;

        virtual void pr_tree(int_type level) = 0;
        virtual bool check_invariant() = 0;
        



};