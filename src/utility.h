#pragma once
#include <cmath>

//debug
#define _ << " " <<
#define NL "\n"
#define V(x) std::string(#x "=") << (x) << " " //"x=...

namespace utility {
    /* easy way to print bit represenation */ 
    void pr_bit(uint64_t x) {
        std::bitset<64> bs(x);
        std::cout << bs << "\n";
    }

    /* vector<bool> -> vector<uint64_t>  */ 
    std::vector<uint64_t> convert_bitstring(std::vector<bool> &bits) {
        uint64_t n = (bits.size() + 63) / 64;
        uint32_t i = 0, j = 0;
        std::vector<uint64_t> v(n);
        for(bool b : bits) {
            v[j] |= (1LL << i) * b;
            i++;
            if(i == 64) {
                i = 0;
                j++;
            }
        }
        return v;
    }
    


    /* To handle base case of below recursive Variadic function Template */ 
    template <uint32_t N = 1>
    void print() {std::cout << "\n";}
    
    /* Variadic function Template that takes variable number of arguments and prints all of them. */
    template <uint32_t N = 1, typename T, typename... Types>
    void print(T var1, Types... var2) {
        std::cout << var1 << std::string(N, ' ');
        print<N>(var2...);
    }

    /* Round to n digits. */
    double round(double x, uint32_t n) {
        double pow_10 = std::pow(10, n);
        return std::round(x * pow_10) / pow_10;
    }
}