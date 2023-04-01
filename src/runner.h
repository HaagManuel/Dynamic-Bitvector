#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "src/dynamic_bv.h"

//debug
#define _ << " " <<
#define NL "\n"
#define V(x) std::string(#x "=") << (x) << " " //"x=...


#define DBG if constexpr(DEBUG_OUT) std::cout
#define ONLY_DBG if constexpr(DEBUG_OUT)

class RunnerBV {    
    private:
        static const bool DEBUG_OUT = false; 
    public:
        RunnerBV(std::string &in, std::string &out) : in_file(in), out_file(out) {};
        void read_file() {
            DBG << "-> reading file from: " << in_file << "\n";
            std::ifstream file;
            file.open(in_file);
            output_len = 0;
            uint32_t n, i, b;
            std::string cmd;
            file >> n;
            while(n--) { file >> b; init_bits.push_back(b); }
            while(file >> cmd) {
                char c = cmd[0]; 
                cmds.push_back(c);
                if(c == 'd' || c == 'f') { //delete, flip
                    file >> i;
                    args.push_back(i);
                }
                else { //insert, rank, select
                    file >> i >> b;
                    args.push_back(i);
                    args.push_back(b);
                }
                output_len += (c == 'r' || c == 's');
            }
            file.close();
            output.resize(output_len);
            DBG << "-> read file \n";
            DBG << "-> init_bits=" << init_bits.size() _ ", num_operations=" << cmds.size() _ ", output_len=" << output_len << "\n";
        }

        template <class BV>
        void run_operations() {
            std::vector<uint64_t> compressed_bits = utility::convert_bitstring(init_bits);
            uint32_t i = 0, j = 0;

            DBG << "-> running operations \n";

            auto t1 = std::chrono::high_resolution_clock::now();
            
            /* construction */
            BV bv(compressed_bits, init_bits.size());

            /* operations */
            for(char c : cmds) {
                switch(c) {
                    case 'i':
                        bv.insert(args[i], args[i + 1]); i += 2;
                        break;
                    case 'd':
                        bv.remove(args[i++]);
                        break;
                    case 'f':
                        bv.flip(args[i++]);
                        break;
                    case 'r':
                        output[j++] = bv.rank(args[i + 1], args[i]); i += 2; // input: [0|1] i, rank(int, bool)
                        break;
                    case 's':
                        output[j++] = bv.select(args[i + 1], args[i]); i += 2; // input: [0|1] i, selct(int, bool)
                        break;
                }
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;

            auto memory_data = bv.memory_consumption();
            uint64_t total_memory_bits = memory_data.memory_total * 8; //bytes -> bits
            ONLY_DBG {memory_data.print();}

            //final output
            std::cout << "RESULT algo=bv name<manuel_haag> time=<" << time << "> space=<" << total_memory_bits << "> \n";
        }

        void output_file() {
            DBG << "-> writing output to file: " << out_file << "\n \n";
            std::ofstream file;
            file.open(out_file);
            for(auto x : output) file << x << "\n";
            file.close();
            return;
        }
    private:
        std::string &in_file;
        std::string &out_file;
        std::vector<char> cmds;
        std::vector<bool> init_bits;
        std::vector<uint32_t> args;
        std::vector<uint32_t> output;
        // MemoryData memory_data;
        // MemoryData memory;
        // DynamicBV::MemoryData mem;
        uint32_t output_len;
};



class RunnerBP {    
    private:
        static const bool DEBUG_OUT = false; 
    public:
        RunnerBP(std::string &in, std::string &out) : in_file(in), out_file(out) {};
        void read_file() {
            DBG << "-> reading file from: " << in_file << "\n";
            std::ifstream file;
            file.open(in_file);
            output_len = 0;
            uint32_t i, v, k;
            std::string cmd;
            while(file >> cmd) {
                char c = cmd[0]; 
                cmds.push_back(c);
                switch(c) {
                    case 'd': //deletenode
                        file >> v;
                        args.push_back(v);
                        break;
                    case 'i': //insertchild
                        file >> v >> i >> k;
                        args.push_back(v);
                        args.push_back(i);
                        args.push_back(k);
                        break;
                    case 'c': //child
                        file >> v >> i;
                        args.push_back(v);
                        args.push_back(i);
                        break;
                    case 's': //subtree_size
                        file >> v;
                        args.push_back(v);
                        break;
                    case 'p': //parent
                        file >> v;
                        args.push_back(v);
                        break;
                }
                output_len += (c == 'r' || c == 's' || c == 'p');
            }
            file.close();
            output.resize(output_len);
            DBG << "-> read file \n";
            DBG << "-> num_operations=" << cmds.size() _ ", output_len=" << output_len << "\n";
        }

        template <class BP>
        void run_operations() {
            uint32_t i = 0, j = 0;
            uint32_t v;
            DBG << "-> running operations \n";

            auto t1 = std::chrono::high_resolution_clock::now();
            
            /* construction */
            BP bp;
            bp.only_root();
            
            /* operations */
            for(char c : cmds) {
                v = bp.select(args[i], true); // v starts at 1
                // std::cout << "cmd " << c  << " " << V(v) _ V(args[i]) _ V(i) _ V(j) "\n";
                switch(c) {
                    case 'd': //deletenode
                        bp.delete_node(v);
                        i += 1;
                        break;
                    case 'i': //insertchild
                        bp.insert_child(v, args[i + 1], args[i + 2]);
                        i += 3;
                        break;
                    case 'c': //child
                        output[j++] = bp.child_i(v, args[i + 1]);
                        i += 2;
                        break;
                    case 's': //subtree_size
                        output[j++] = bp.subtree_size(v);
                        i += 1;
                        break;
                    case 'p': //parent
                        output[j++] = bp.parent(v);
                        i += 1;
                        break;
                }
            }

            //computing out degree
            std::vector<uint32_t> out_deg = bp.all_out_deg();
            for(auto x : out_deg) output.push_back(x);
            
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;


            auto memory_data = bp.memory_consumption();
            uint64_t total_memory_bits = memory_data.memory_total * 8; //bytes -> bits
            ONLY_DBG {memory_data.print();}

            //final output
            std::cout << "RESULT algo=bp name<manuel_haag> time=<" << time << "> space=<" << total_memory_bits << "> \n";
        }

        void output_file() {
            DBG << "-> writing output to file: " << out_file << "\n \n";
            std::ofstream file;
            file.open(out_file);
            for(auto x : output) file << x << "\n";
            file.close();
            return;
        }
    private:
        std::string &in_file;
        std::string &out_file;
        std::vector<char> cmds;
        std::vector<bool> init_bits;
        std::vector<uint32_t> args;
        std::vector<uint32_t> output;
        // MemoryData memory_data;
        // MemoryData memory;
        // DynamicBV::MemoryData mem;
        uint32_t output_len;
};