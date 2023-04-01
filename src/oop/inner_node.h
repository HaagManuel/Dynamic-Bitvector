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

class InnerNode final : public Node{
    public:
        using int_type = uint32_t;

        InnerNode(Node* l, Node* r) : left(l), right(r), num(l -> get_num()), ones(l -> get_ones()), height(1) {};
        ~InnerNode() {
            left -> ~Node();
            right -> ~Node();
        }

        bool access(int_type i) override {
            if(i < num) { return left -> access(i);}
            else         { return right -> access(i - num);}
        }
        bool flip(int_type i) override {
            bool b;
            if(i < num) {
                b = left -> flip(i);
                ones += b? 1 : -1;
            }
            else { b = right -> flip(i - num);}
            return b;
        }
        Node* insert(int_type i, bool b) override {
            if(i <= num) {
                num++; ones += b;
                left = left -> insert(i, b);
            }
            else { right = right -> insert(i - num, b);}
            Node* to_return = this;
            if(is_imbalanced()) to_return = rebalance();
            update_height();
            return to_return;
        }

        int_type rank1(int_type i, int_type ones_before) override {
            if(i <= num) { return left -> rank1(i, ones_before); }
            else         { return right -> rank1(i - num, ones_before + ones);}
        }
        int_type select1(int_type i, int_type bits_before) override {
            if(i <= ones) { return left -> select1(i, bits_before);}
            else          { return right -> select1(i - ones, bits_before + num);}
        }
        int_type select0(int_type i, int_type bits_before) override {
            int_type zeros = num - ones;
            if(i <= zeros) { return left -> select0(i, bits_before);}
            else          {  return right -> select0(i - zeros, bits_before + num);}
        }

        int_type get_num() override { return num;} 
        int_type get_ones() override { return ones;} 
        int_type get_height() override { return height;}

        void pr_tree(int_type level) override {
            std::cout << std::string(level, ' ' ) << "( height, num, ones:" 
            << " " << height 
            << " " << num 
            << " " << ones 
            << "\n"; 
            left -> pr_tree(level + 1);
            right -> pr_tree(level + 1);
            std::cout << std::string(level, ' ' ) << ")" << "\n";
        } 

        bool check_invariant() override {
            int32_t balance = get_balance();
            if(std::abs(balance) > 1) return false;
            if(!left -> check_invariant()) return false;
            if(!right -> check_invariant()) return false;
            return true;
        }

        int32_t get_balance() override {
            return left -> get_height() - right -> get_height();
        }

        bool is_imbalanced() {
            return std::abs(get_balance()) > 1;
        }

        void update_height() override {
            height = 1 + std::max(left -> get_height(), right -> get_height());
        }
        
        void set_num(int_type n) override {  
            num = n;
        }
        
        void set_ones(int_type o) override {  
            ones = o;
        }
        
        /*
             y                               x
            / \     Right Rotation          /  \
           x   T3   - - - - - - - >        T1   y 
          / \       < - - - - - - -            / \
         T1  T2     Left Rotation            T2  T3
        */
        Node* right_rotate(Node* z) {
            InnerNode* y = (InnerNode*) z;
            InnerNode* x = (InnerNode*) y -> left; //must be inner node by invariant
            Node *T2 = x -> right;
            
            //rotation
            x -> right = y;
            y -> left = T2;

            //update ones
            y -> set_ones(y -> get_ones() - x -> get_ones());
            y -> set_num(y -> get_num() - x -> get_num());

            y -> update_height();
            x -> update_height();
            return x; //new root
        }

        Node* left_rotate(Node* z) {
            InnerNode* x = (InnerNode*) z;
            InnerNode* y = (InnerNode*) x -> right; //must be inner node by invariant
            Node* T2 = y -> left;
            
            //rotation
            x -> right = T2;
            y -> left = x;

            //update ones
            y -> set_ones(y -> get_ones() + x -> get_ones());
            y -> set_num(y -> get_num() + x -> get_num());

            x -> update_height();
            y -> update_height();
            return y; //new root
        }

        //returns new subtree root
        Node* rebalance() { 
            assert(is_imbalanced()); //otherwise returns NULL
            int32_t balance = get_balance();
            Node *subtree_root = NULL;
            if(balance > 1) { //left
                int32_t balance_y = left -> get_balance();
                if(balance_y > 0) { //left left
                    subtree_root = right_rotate(this);
                }
                else { //left right
                    left = left_rotate(left);
                    subtree_root = right_rotate(this);
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
                int32_t balance_y = right -> get_balance();    
                if(balance_y < 0) { //right right           
                    subtree_root = left_rotate(this);
                }
                else {
                    right = right_rotate(right);
                    subtree_root = left_rotate(this);
                }
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

    private:
        Node* left;
        Node* right;
        int_type num;
        int_type ones;
        int_type height;

};