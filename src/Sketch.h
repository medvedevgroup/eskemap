#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <cstdint>
#include <utility>
#include <list>
#include <string>
#include <iostream>

//The k-mer size
#define K 9
//The hash value threshold
#define MAX_HASH 26214

//A Sketch is a list of offset hash pairs
using Sketch = std::list<std::pair<uint32_t, uint64_t>>;

//This function returns a numerical value for each nucleotide base
inline uint64_t getBaseNb(const char& c){
	switch(c){
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			std::cerr << "WARNING: Unsupported character in input sequence detected!" << std::endl;
			return 0;
	}
}

//This function calculates a hash for the given numerical representation of a k-mer and a mask of the form (4**k)-1 where k is the k-mer length;
//it originates from the minimap2 source code (the function hash64 in sketch.c)
static inline uint64_t getHash(uint64_t kmer, uint64_t mask){
	//Testing
	std::cout << "~kmer: " << ~kmer << std::endl;
	
	kmer = (~kmer + (kmer << 21)) & mask; // kmer = (kmer << 21) - kmer - 1;

	//Testing
	std::cout << kmer << std::endl;

	kmer = kmer ^ kmer >> 24;
	kmer = ((kmer + (kmer << 3)) + (kmer << 8)) & mask; // kmer * 265
	kmer = kmer ^ kmer >> 14;
	kmer = ((kmer + (kmer << 2)) + (kmer << 4)) & mask; // kmer * 21
	kmer = kmer ^ kmer >> 28;
	kmer = (kmer + (kmer << 31)) & mask;
	return kmer;
}

//This function calculates the numerical representation of a k-mer
uint64_t calcKmerNb(const std::string& kmer);

//This function builds the sketch of a sequence using a given hash threshold
Sketch buildSketch(const std::string& s, const uint64_t& ht);

#endif