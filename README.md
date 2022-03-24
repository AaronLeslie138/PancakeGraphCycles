Cycles.cpp finds the lengths at which cycles appear in pancake graphs.

Core algorithm is modification of hawick's circuit and loops enumeration algorithm, as implemented in the boost library, to be multi-threaded and to self-impose depth-limits when cycles have been found for given depths (as this project is only interested in confirming the lengths at which cycles appear, not enumerating them.)