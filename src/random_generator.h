#pragma once
#include <random>
#include <algorithm>

namespace random_gen {
    size_t random_index(const size_t upper) {
        static thread_local std::mt19937 generator;
        std::uniform_int_distribution<int> distribution(0, upper - 1);
        return distribution(generator);
    }
    //inclusive in both directions
    size_t random_range(const size_t lower, const size_t upper) {
        static thread_local std::mt19937 generator;
        std::uniform_int_distribution<int> distribution(lower, upper);
        return distribution(generator);
    }

    bool random_bool() {
        return random_index(2) == 0;
    }

    // Prob = 1/p, p in [0,1]
    bool prob_p(const double &p) {
        static thread_local std::minstd_rand gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(0, 1);
        return dist(gen) <= p;
    }

    std::vector<bool> random_bit_string(uint64_t n) {
        std::vector<bool> v(n);
        for(uint64_t i = 0; i < n; i++) {
            v[i] = random_bool();
        }
        return v;
    }

    template <class T>
    void shuffle(std::vector<T> &v) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(v.begin(), v.end(), g);
    }
    
    std::vector<uint32_t> random_permutation(uint32_t n) {
        std::vector<uint32_t> v(n);
        std::iota(v.begin(), v.end(), 0);
        shuffle<uint32_t>(v);
        return v;
    }

    
}