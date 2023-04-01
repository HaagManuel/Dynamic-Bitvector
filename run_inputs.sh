cd release
# cd debug
make
./main bv ../inputs/bsp_bv/example_bv_10k.txt ../outputs/bv_out-example_bv_10k.txt
# ./main bv ../inputs/bsp_bv/example_bv_100k.txt ../outputs/bv_out-example_bv_100k.txt #-> input contains artefact!
./main bv ../inputs/bsp_bv/example_bv_1M.txt ../outputs/bv_out-example_bv_1M.txt

./main bp ../inputs/bsp_tree/example_tree_d6_c5-10.txt ../outputs/bp_out-example_tree_d6_c5-10.txt
./main bp ../inputs/bsp_tree/example_tree_d8_c5.txt ../outputs/bp_out-example_tree_d8_c5.txt
./main bp ../inputs/bsp_tree/example_tree_d9_c4.txt ../outputs/bp_out-example_tree_d9_c4.txt


diff ../inputs/bsp_bv/example_bv_10k_output.txt ../outputs/bv_out-example_bv_10k.txt
# diff ../inputs/bsp_bv/example_bv_100k_output.txt ../outputs/bv_out-example_bv_100k.txt
diff ../inputs/bsp_bv/example_bv_1M_output.txt ../outputs/bv_out-example_bv_1M.txt

diff ../inputs/bsp_tree/example_tree_d6_c5-10_output.txt ../outputs/bp_out-example_tree_d6_c5-10.txt
diff ../inputs/bsp_tree/example_tree_d8_c5_output.txt ../outputs/bp_out-example_tree_d8_c5.txt
diff ../inputs/bsp_tree/example_tree_d9_c4_output.txt ../outputs/bp_out-example_tree_d9_c4.txt



