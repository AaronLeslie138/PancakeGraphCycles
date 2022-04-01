Cycles.cpp finds the lengths at which cycles appear in pancake graphs.

Core algorithm is modification of hawick's circuit and loops enumeration algorithm, as implemented in the boost library, to be multi-threaded and to self-impose depth-limits when cycles have been found for given depths (as this project is only interested in confirming the lengths at which cycles appear, not enumerating them.)

## How to Run
Download Boost v1.73.0 to project directory, without replacing modified hawick_circuits.hpp supplied in project.
Compile with GCC
Run program and supply sides and permutations
Will generate folder corresponding to specific graph, with files containing a sample cycle for each length at which a cycle was found

##Directed Graphs
By default, computes for undirected pancake graphs. 
Remove undirectedS modified from adjancy list on line 170 to run on directed graphs
> typedef boost::adjacency_list< boost::vecS, boost::vecS/*, boost::undirectedS*/ > Graph;
(May also want to replace instances of "-U" with "-D" in main, as "-U" is intended to notate undirect graphs.
