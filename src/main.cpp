#include <iostream>
#include <cassert>
#include <array>
#include <stack>
#include <cmath>
#include <iomanip>   

#include <chrono>
#include <variant>
#include <vector>

#include "src/random_generator.h"
#include "src/naive_bv.h"
#include "src/static_small_bv.h"
#include "src/dynamic_bv.h"
#include "src/oop/dynamic_bv_oop.h"
#include "src/variant/dynamic_bv_variant.h"
#include "src/runner.h"
#include "src/utility.h"

#include "src/dynamic_tree_bp.h"

#include "src/dynamic_small_bv.h"

template<class BV> 
struct RandomOperation {
    RandomOperation(BV &bv) : _bv(bv), cur_size(0) {
        count[0] = 0; count[1] = 0;
    };

    //n = num total repetitions, m = times insert and remove operation is perfomed 

    void random_insert(uint32_t n = 1, uint32_t m = 1) {
        uint64_t x = 0;
        for(uint32_t i = 0; i < n; i++) {
            uint32_t j = random_gen::random_index(cur_size + 1);   
            bool b = random_gen::random_bool();
            for(uint32_t k = 0; k < m; k++) {
                _bv.insert(j, b);
                x += _bv.access(j);
                x &= ((1LL << 60) - 1);
                count[b]++;
                cur_size++;
            }
        } 
        if(x == 12304123) std::cout << "hurra! \n"; //force evaluation
    }
    void random_remove(uint32_t n = 1, uint32_t m = 1) {
        for(uint32_t i = 0; i < n; i++) {
            uint32_t j = random_gen::random_index(cur_size);   
            for(uint32_t k = 0; k < m; k++) {
                bool bit = _bv.access(j);
                count[0] -= bit; count[1] -= bit;
                cur_size--;
                _bv.remove(j);
                if(j == cur_size) j--;
                // std::cout << V(j) _ V(cur_size) << "\n";
            }
        }
    }
    void random_flip(uint32_t n = 1) {
        for(uint32_t i = 0; i < n; i++) {
            uint32_t j = random_gen::random_index(cur_size);   
            _bv.flip(j);
            int32_t diff = _bv.access(j) ? 1 : -1;
            count[0] -= diff; count[1] += diff;
        }
    }

    void random_access_rank(uint32_t n = 1) {
        uint64_t x = 0;
        for(uint32_t i = 0; i < n; i++) {
            uint32_t j = random_gen::random_index(cur_size);   
            x += _bv.access(j);
            x += _bv.rank(j, (j % 2) == 0);
            x &= ((1LL << 60) - 1);
        }
        if(x == 12304123) std::cout << "hurra! \n"; //force evalutaiton
    }

    void all_access_rank() {
        uint64_t x = 0;
        for(uint32_t i = 0; i < cur_size; i++) {
            x += _bv.access(i);
            x += _bv.rank(i, true);
            x += _bv.rank(i, false);
            x &= ((1LL << 60) - 1);
        }
        if(x == 12304123) std::cout << "hurra! \n"; //force evalutaiton
    }

    void permuted_all_access_rank() {
        uint64_t x = 0;
        std::vector<uint32_t> v = random_gen::random_permutation(cur_size);
        for(uint32_t i = 0; i < cur_size; i++) {
            uint32_t j = v[i];
            x += _bv.access(j);
            x += _bv.rank(j, true);
            x += _bv.rank(j, false);
            x &= ((1LL << 60) - 1);
        }
        if(x == 12304123) std::cout << "hurra! \n"; //force evalutaiton
    }

    void random_select(uint32_t n = 1) {
        uint64_t x = 0;
        for(uint32_t i = 0; i < n; i++) {
            uint32_t j0 = 1 + random_gen::random_index(count[i % 2]);   
            x += _bv.access(j0);
            x += _bv.select(j0, i % 2);
            x &= ((1LL << 60) - 1);
        }
        if(x == 12304123) std::cout << "hurra! \n"; //force evalutaiton
    }
    BV &_bv;
    std::vector<uint32_t> indices;
    std::vector<bool> bools;
    uint32_t cur_size;
    uint32_t count[2];
};


template<class BV> 
std::tuple<double, double, double, double> benchmark(uint32_t N, uint32_t runs) {
    double avg_insert = 0, avg_query = 0, avg_flip = 0, avg_remove = 0;
    for(uint32_t i = 0; i < runs; i++) {
        BV bv; RandomOperation op(bv);
        /* insert */
        auto t1 = std::chrono::high_resolution_clock::now();
        op.random_insert(N);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        avg_insert += time;

        /* flip */
        t1 = std::chrono::high_resolution_clock::now();
        // op.random_flip(N);
        t2 = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        avg_flip += time;

        /* query */
        t1 = std::chrono::high_resolution_clock::now();
        // op.random_access_rank(N);
        // op.random_select(N);
        t2 = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        avg_query += time;
        
        // if(i == 0) {
        //     auto memory =  bv.memory_consumption();
        //     memory.print();
        // }
        /* remove */
        t1 = std::chrono::high_resolution_clock::now();
        op.random_remove(N);
        t2 = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        avg_remove += time;

    }
    avg_insert /= (double) runs;
    avg_flip /= (double) runs;
    avg_query /= (double) runs;
    avg_remove /= (double) runs;

    uint32_t n = 2;
    return {utility::round(avg_insert, n), utility::round(avg_flip, n), utility::round(avg_query, n), utility::round(avg_remove, n)};
}

template<class BV> 
void benchmark2 (std::vector<uint64_t> stops, uint32_t N, uint32_t M, uint32_t runs) {
    utility::print<10>("stop", "insert", "flip", "rank", "select", "remove");
    BV bv; RandomOperation op(bv);
    for(auto stop : stops) {
        double avg_insert = 0, avg_rank = 0, avg_select = 0, avg_flip = 0, avg_remove = 0;
        uint64_t diff = stop - bv.size();
        op.random_insert(diff, 1); //bring to bv to desired size

        for(uint32_t i = 0; i < runs; i++) {
            /* insert */
            auto t1 = std::chrono::high_resolution_clock::now();
            op.random_insert(N, M);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            avg_insert += time;

            /* flip */
            t1 = std::chrono::high_resolution_clock::now();
            op.random_flip(N * M);
            t2 = std::chrono::high_resolution_clock::now();
            time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            avg_flip += time;

            /* rank */
            t1 = std::chrono::high_resolution_clock::now();
            op.random_access_rank(N * M);
            t2 = std::chrono::high_resolution_clock::now();
            time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            avg_rank += time;
            
            /* select */
            t1 = std::chrono::high_resolution_clock::now();
            op.random_select(N * M);
            t2 = std::chrono::high_resolution_clock::now();
            time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            avg_select += time;
            
            /* remove */
            t1 = std::chrono::high_resolution_clock::now();
            op.random_remove(N, M);
            t2 = std::chrono::high_resolution_clock::now();
            time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            avg_remove += time;

        }
        avg_insert /= (double) runs;
        avg_flip /= (double) runs;
        avg_rank /= (double) runs;
        avg_select /= (double) runs;
        avg_remove /= (double) runs;

        uint32_t n = 2;
        avg_insert = utility::round(avg_insert, n);
        avg_flip = utility::round(avg_flip, n);
        avg_rank = utility::round(avg_rank, n);
        avg_select = utility::round(avg_select, n);
        avg_remove = utility::round(avg_remove, n);
        utility::print<10>(stop, avg_insert, avg_flip, avg_rank, avg_select, avg_remove);
    }
    
    return;
}


template<class BV> 
void benchmark_build_up_only_insert (std::vector<uint64_t> stops, uint32_t runs) {
    utility::print<10>("#algo", "run", "n", "insert");
    for(uint32_t j = 0; j < runs; j++) {
        BV bv; RandomOperation op(bv);
        std::vector<double> time_insert(stops.size(), 0);
        uint32_t k = 0;
        for(auto stop : stops) {
            uint64_t diff = stop - bv.size();
            auto t1 = std::chrono::high_resolution_clock::now();
            op.random_insert(diff, 1);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            assert(bv.size() == stop);
            time_insert[k] += time;
            k++;
        }
        for(uint32_t i = 0; i < stops.size(); i++) {
            utility::print<10>("bv_oop", j, stops[i], time_insert[i]);
        }
    }   
    return;
}


template<class BV> 
void benchmark_build_up (std::vector<uint64_t> stops, uint32_t runs, bool header = true) {
    if(header) utility::print<10>("#algo", "run", "n", "leaf_blocks", "insert", "remove", "memory_leafs", "memory_nodes", "memory_total", "fraction_stored_used_bits", "bits_per_stored_bit");
    for(uint32_t j = 0; j < runs; j++) {
        BV bv; RandomOperation op(bv);
        std::vector<double> time_insert(stops.size(), 0), time_remove(stops.size(), 0);
        std::vector<typename BV::MemoryData> memorydata;
        uint32_t k = 0;
        for(auto stop : stops) {
            uint64_t diff = stop - bv.size();
            auto t1 = std::chrono::high_resolution_clock::now();
            op.random_insert(diff, 1);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            assert(bv.size() == stop);
            time_insert[k] += time;
            memorydata.push_back(bv.memory_consumption());
            k++;
        }
        k = 0;
        uint64_t removed = 0;
        for(auto stop : stops) {
            // std::cout << "rmv bv " << V(stop) << "\n";
            uint64_t diff = stop - removed;
            removed += diff;

            auto t1 = std::chrono::high_resolution_clock::now();
            op.random_remove(diff, 1);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            
            time_remove[k] += time;
            k++;
        }
        assert(bv.size() == 0);

        for(uint32_t i = 0; i < stops.size(); i++) {
            auto mem = memorydata[i];
            BV bv;
            utility::print<10>("bv", j, stops[i], bv.leaf_blocks(),  time_insert[i], time_remove[i],
             (double)mem.memory_leafs, (double)mem.memory_nodes, (double)mem.memory_total, mem.fraction_stored_used_bits, mem.bits_per_stored_bit
            );
        }
    }   
    return;
}



//insert = insert_child (leaf), remove = delete_node
template<class BP> 
void benchmark_build_up_bp (std::vector<uint64_t> stops, uint32_t runs, bool header = true) {
    if(header) utility::print<9>("#algo", "run", "n", "leaf_blocks", "insert", "remove", "memory_leafs", "memory_nodes", "memory_total", "fraction_stored_used_bits", "bits_per_stored_bit");
    for(uint32_t j = 0; j < runs; j++) {
        BP bp; bp.only_root();
        std::vector<double> time_insert(stops.size(), 0), time_remove(stops.size(), 0);
        std::vector<typename BP::MemoryData> memorydata;
        uint32_t k = 0;
        for(auto stop : stops) {
            uint64_t diff = stop / 2 - bp.size() / 2;
            auto t1 = std::chrono::high_resolution_clock::now();
            //insert random leafs
            for(uint i = 0; i < diff; i++) {
                int32_t k = 1 + random_gen::random_index(i + 1); //random node, node number start at 0
                int32_t v = bp.select(k, true);
                bp.insert_child(v, 1, 0);
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            assert(bp.size() == stop);
            time_insert[k] += time;
            memorydata.push_back(bp.memory_consumption());
            k++;
        }
        //print size
        
        k = 0;
        uint64_t removed = 0;
        stops.back() = stops.back() - 2; //avoid removing root
        for(auto stop : stops) {
            uint64_t diff = stop / 2 - removed;
            removed += diff;

            auto t1 = std::chrono::high_resolution_clock::now();
            //remove random nodes, but not root
            for(uint i = 0; i < diff; i++) {
                int32_t k = random_gen::random_range(2, bp.size() / 2);
                int32_t v = bp.select(k, true);
                bp.delete_node(v);
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            time_remove[k] += time;
            k++;
        }
        assert(bp.size() == 2);

        for(uint32_t i = 0; i < stops.size(); i++) {
            auto mem = memorydata[i];
            BP bp;
            utility::print<10>("bp", j, stops[i], bp.leaf_blocks(),  time_insert[i], time_remove[i],
             (double)mem.memory_leafs, (double)mem.memory_nodes, (double)mem.memory_total, mem.fraction_stored_used_bits, mem.bits_per_stored_bit
            );
        }
    }   
    return;
}



//builds up data structure, at every stop num_flips flips are performed and measured
template<class BV> 
void benchmark_flip (std::vector<uint64_t> stops, uint32_t runs, uint64_t num_flips, std::string algo_name, bool header = true) {
    if(header) utility::print<8>("#algo", "run", "n", "leaf_blocks", "flip");
    for(uint32_t j = 0; j < runs; j++) {
        BV bv; RandomOperation op(bv);
        std::vector<double> time_flip(stops.size(), 0);
        uint32_t k = 0;
        for(auto stop : stops) {
            uint64_t diff = stop - bv.size();
            op.random_insert(diff, 1);
            auto t1 = std::chrono::high_resolution_clock::now();
            op.random_flip(num_flips);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
            assert(bv.size() == stop);
            time_flip[k] += time;
            k++;
        }
        for(uint32_t i = 0; i < stops.size(); i++) {
            BV bv;
            utility::print<10>(algo_name, j, stops[i], bv.leaf_blocks(), time_flip[i]);
        }
    }   
    return;
}



//n = num total repetitions, m = times insert and remove operation is perfomed 

void test() {
    //test in release mode
    // uint32_t N = 1e7, runs = 5;
    // uint32_t N = 1e7, runs = 1;
    uint32_t N = 1e5, runs = 1;
    double insert, query, flip, remove;
    auto pr = [&]() {
        utility::print<10>(insert, flip, query, remove);
    };

    std::cout << V(N) _ V(runs) << "\n";
    std::cout << "BV \n";
    utility::print<5>("insert [ms]", "flip [ms]", "query [ms]", "remove [ms]");
    // std::cout << "insert [ms]    flip [ms]    query [ms]    remove [ms]   " << V(N) << "    " << V(runs) << "\n";
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBV<8>>(N, runs);
    // pr();
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBV<2>>(N, runs);
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBV<4>>(N, runs);
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBV<16>>(N, runs);
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBV<32>>(N, runs);
    std::tie(insert, flip, query, remove) = benchmark<DynamicBV<64>>(N, runs);
    pr();
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBV<128>>(N, runs);
    // pr();


    // std::cout << "BV OOP \n";
    // std::cout << "insert [ms]    flip [ms]    query [ms]    " << V(N) << "    " << V(runs) << "\n";
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBVOop<4>>(N, runs);
    // std::cout << insert  << "    " << flip << "    " << query << "\n";
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBVOop<8>>(N, runs);
    // std::cout << insert  << "    " << flip << "    " << query << "\n";
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBVOop<32>>(N, runs);
    // std::cout << insert  << "    " << flip << "    " << query << "\n";

    // std::cout << "BV Variant \n";
    // std::cout << "insert [ms]    flip [ms]    query [ms]    " << V(N) << "    " << V(runs) << "\n";
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBVVariant<4>>(N, runs);
    // std::cout << insert  << "    " << flip << "    " << query << "\n";
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBVVariant<8>>(N, runs);
    // std::cout << insert  << "    " << flip << "    " << query << "\n";
    // std::tie(insert, flip, query, remove) = benchmark<DynamicBVVariant<32>>(N, runs);
    // std::cout << insert  << "    " << flip << "    " << query << "\n";
}

void test2() {
    uint64_t N = 10000, M = 1, last_break = 1e6, runs = 5;
    std::vector<uint64_t> stops;
    for(uint64_t i = 1; i < last_break; i *= 10) stops.push_back(i);
    benchmark2<DynamicBV<64>>(stops, N, M, runs);


}

//1 run, 2^28, ~ 5min
//1 run, 2^27, ~ 2min
//1 run, 2^26, ~ 1min 6:27

//bp 1 run 2^20 -> 16000 ms insert, 
void test_build_up() {
    uint64_t last_break = (1 << 20), runs = 5;
    std::vector<uint64_t> stops;
    for(uint64_t i = 2; i <= last_break; i *= 2) stops.push_back(i);
    // benchmark_build_up<DynamicBV<64>>(stops, runs, true);
    benchmark_build_up_bp<DynamicTreeBP<64>>(stops, runs, false);
}

void test_build_up_by_leafs() {
    uint64_t last_break = (1 << 18), runs = 5;
    std::vector<uint64_t> stops = {last_break};
    benchmark_build_up<DynamicBV<4>>(stops, runs, true);
    benchmark_build_up<DynamicBV<8>>(stops, runs, false);
    benchmark_build_up<DynamicBV<16>>(stops, runs, false);
    benchmark_build_up<DynamicBV<32>>(stops, runs, false);
    benchmark_build_up<DynamicBV<64>>(stops, runs, false);
    benchmark_build_up<DynamicBV<128>>(stops, runs, false);
    benchmark_build_up<DynamicBV<256>>(stops, runs, false);

    benchmark_build_up_bp<DynamicTreeBP<4>>(stops, runs, false);
    benchmark_build_up_bp<DynamicTreeBP<8>>(stops, runs, false);
    benchmark_build_up_bp<DynamicTreeBP<16>>(stops, runs, false);
    benchmark_build_up_bp<DynamicTreeBP<32>>(stops, runs, false);
    benchmark_build_up_bp<DynamicTreeBP<64>>(stops, runs, false);
    benchmark_build_up_bp<DynamicTreeBP<128>>(stops, runs, false);
    benchmark_build_up_bp<DynamicTreeBP<256>>(stops, runs, false);

}

void test_flips() {
    uint64_t last_break = (1 << 20), runs = 5, num_flips = 1e6;
    std::vector<uint64_t> stops;
    for(uint64_t i = 2; i <= last_break; i *= 2) stops.push_back(i);
    // benchmark_flip<DynamicBV<64>>(stops, runs, num_flips, "flip_it_opt", true);
    // benchmark_flip<DynamicBV<64>>(stops, runs, num_flips, "flip_it", false);
    // benchmark_flip<DynamicBV<64>>(stops, runs, num_flips, "flip_rec", false);
    benchmark_flip<DynamicBVOop<64>>(stops, runs, num_flips, "flip_rec_oop", false);
}


template<class BV> 
void construction_benchmark(uint64_t N, uint64_t runs) {
    double naive_front = 0, naive_back = 0, fast = 0, fast_with_compress = 0;;
    for(uint32_t i = 0; i < runs; i++) {
        std::vector<bool> seq = random_gen::random_bit_string(N);
        std::vector<uint64_t> compressed = utility::convert_bitstring(seq);

        /* front */
        auto t1 = std::chrono::high_resolution_clock::now();
        BV bv_naive_front;
        for(bool b : seq) bv_naive_front.insert(0, b);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        naive_front = time;

        /* back */
        t1 = std::chrono::high_resolution_clock::now();
        BV bv_naive_back;
        for(uint64_t i = 0; i < N; i++) bv_naive_front.insert(i, seq[i]);
        t2 = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        naive_back = time;

        /* fast */
        t1 = std::chrono::high_resolution_clock::now();
        BV bv_fast(compressed, N);
        t2 = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        fast = time;

        /* fast with compression time*/
        t1 = std::chrono::high_resolution_clock::now();
        std::vector<uint64_t> compressed2 = utility::convert_bitstring(seq);
        BV bv_fast2(compressed2, N);
        t2 = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000.;
        fast_with_compress = time;

        BV bv;
        utility::print<8>("bv", i, N, bv.leaf_blocks(), naive_front, naive_back, fast, fast_with_compress);
    }
    
    
}

void test_construction() {
    uint64_t last_break = (1 << 20), runs = 5;
    std::vector<uint64_t> stops;
    for(uint64_t i = 2; i <= last_break; i *= 2) stops.push_back(i);

    utility::print<6>("#run", "n", "leaf_blocks", "front", "back", "rec", "compress_&_rec");
    for(auto x : stops) {
        construction_benchmark<DynamicBV<4>>(x, runs);
        construction_benchmark<DynamicBV<8>>(x, runs);
        construction_benchmark<DynamicBV<16>>(x, runs);
        construction_benchmark<DynamicBV<32>>(x, runs);
        construction_benchmark<DynamicBV<64>>(x, runs);
        construction_benchmark<DynamicBV<128>>(x, runs);
        construction_benchmark<DynamicBV<256>>(x, runs);
    }
}

void run_from_file_bv() {
    std::string in_dir = "../inputs/bsp_bv";
    // std::string in_file = "example_bv_10k.txt";
    // std::string in_file = "example_bv_100k.txt";
    std::string in_file = "example_bv_1M.txt";
    std::string in_path = in_dir +  "/" +  in_file;

    std::string out_dir = "../outputs";
    std::string out_file = "out_" +  in_file;
    std::string out_path = out_dir +  "/" +  out_file;

    RunnerBV runner(in_path, out_path);
    runner.read_file();
    runner.run_operations<DynamicBV<64>>();
    runner.output_file();
}

void run_from_file_bp() {
    std::string in_dir = "../inputs/bsp_tree";
    // std::string in_file = "example_tree_d6_c5-10.txt";
    std::string in_file = "example_tree_d8_c5.txt"; 
    // std::string in_file = "example_tree_d9_c4.txt"; //to long
    std::string in_path = in_dir +  "/" +  in_file;

    std::string out_dir = "../outputs";
    std::string out_file = "out_" +  in_file;
    std::string out_path = out_dir +  "/" +  out_file;

    RunnerBP runner(in_path, out_path);
    runner.read_file();
    runner.run_operations<DynamicTreeBP<64>>(); 
    runner.output_file();
}

void parse_cmdline(int argc, char *argv[]) {
    if(argc != 4) {
        std::cout << "input format: adsprogramma [bv|bp] Eingabedatei Ausgabedatei \n";
        return;
    }
    std::string bv_bp = argv[1];
    std::string file_in = argv[2];
    std::string file_out = argv[3];

    if(bv_bp == std::string("bv")) {
        RunnerBV runner(file_in, file_out);
        runner.read_file();
        runner.run_operations<DynamicBV<64>>();
        runner.output_file();
    }
    else if(bv_bp == std::string("bp")) {
        RunnerBP runner(file_in, file_out);
        runner.read_file();
        runner.run_operations<DynamicTreeBP<64>>();
        runner.output_file();
    }
    else {
        std::cout << "need to specify bv or bp";
    }
}


void vector_allocation_test() {
    int N = 30;
    std::vector<int> v;

    auto print_v = [&]() {
        for(auto x : v) std::cout << x << " ";
        std::cout << "\n";
    };
    for(int i = 0; i < N / 2; i++) {
        v.reserve(v.size() + 2);
        v.resize(v.size() + 2, 33);
        v[v.size() - 1] = i;
        // v[v.size() - 2] = i;
        print_v();
        std::cout << V(v.size()) _ V(v.capacity()) << "\n";
    }
    for(int i = 0; i < N / 2; i++) {
        std::vector<int>(v.begin(), std::max(v.begin() + 1, v.end() - 2)).swap(v); //on element less
        print_v();
        std::cout << V(v.size()) _ V(v.capacity()) << "\n";
    }
}

int main(int argc, char *argv[]) {
    if(argc > 1) {
        parse_cmdline(argc, argv);
        return 0;
    }

    test2();
    test_build_up();
    test_build_up_by_leafs();
    test_flips();
    test_construction();

    vector_allocation_test();

    return 0;
}
