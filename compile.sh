g++ -pthread -w -std=c++11 -O3 -DNDEBUG -DDEGORDER -DFIRSTASC -DSECONDASC  -I /extdata3/parklab/dna/include -L /extdata3/parklab/dna/lib DAF_div.cpp -o daf_10sec -lsdsl -ldivsufsort -ldivsufsort64
