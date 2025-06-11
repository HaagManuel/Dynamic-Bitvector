# Dynamic Bit Vector and Balanced Parenthesis
This repository contains my final project of the [Advanced Datastructure SS2022](https://ae.iti.kit.edu/4264.php) course of ITI Algorithm Engineering at Karlsruher Institute of Technology (KIT).
We implemented a [Dynamic Bit Vector](https://ae.iti.kit.edu/download/kurpicz/2022_advanced_data_structures/03_dynamic_bit_vectors_trees_handout_ss22.pdf) based on balanced search trees (e.g. [AVL-Tree](https://en.wikipedia.org/wiki/AVL_tree)), which supports *access, insert, delete* and *flip* operations of bits.
Additionally, the *rank* operation (number of 0s or 1s before an index $i$) and the *select* operation (position of the $i$-th 0 or 1) are supported.
These operations are crucial when implementing succint, i.e. very memory-efficient, data structures. 

The second task was to implement a succint tree data structure called *Balanced Parentheses* based on *Dynamic Bit Vectors*. 0 and 1 are interpreted as '(' and ')'. 
Parentheses expression are used to represent subtrees in the underlying tree.

**Example of a Dynamic Bit Vector**
![dynamic_bit vector](/images/dynamic_bitvector.png)
[image source](https://ae.iti.kit.edu/download/kurpicz/2022_advanced_data_structures/03_dynamic_bit_vectors_trees_handout_ss22.pdf)



# Compile
```
mkdir release && cd release
cmake ..
```

# Run
```
cd release 
./main [bp|bv] in_file out_file
```



