#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <cstdint>
#include <utility>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "../../software/minimap2/minimap.h"
#include "../../software/minimap2/index.c"

//The k-mer size
#define K 9
//The hash value threshold
#define MAX_HASH 26214
//The FracMinHash ratio
#define HASH_RATIO 0.1
//Size of the considered alphabet
#define ALPHABET_SIZE 4
//Number of times we expect p to occur in t
#define P_MULTIPLICITY 2

using namespace std;
//A PairSketch is a list of offset hash pairs
using PairSketch = list<pair<uint32_t, uint64_t>>;
//A Sketch is a list of hashes
using Sketch = vector<uint64_t>;

//A compare function to sort hashes in a sketch
inline const bool smHshsFrst(const pair<uint32_t, uint64_t>& left, const pair<uint32_t, uint64_t>& right){ return left.second < right.second; }

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
			cerr << "WARNING: Unsupported character in input sequence detected!" << endl;
			return 0;
	}
}

//This function calculates a hash for the given numerical representation of a k-mer and a mask of the form (4**k)-1 where k is the k-mer length;
//it originates from the minimap2 source code (the function hash64 in sketch.c)
static inline uint64_t getHash(uint64_t kmer, uint64_t mask){
	kmer = (~kmer + (kmer << 21)) & mask; // kmer = (kmer << 21) - kmer - 1;
	kmer = kmer ^ kmer >> 24;
	kmer = ((kmer + (kmer << 3)) + (kmer << 8)) & mask; // kmer * 265
	kmer = kmer ^ kmer >> 14;
	kmer = ((kmer + (kmer << 2)) + (kmer << 4)) & mask; // kmer * 21
	kmer = kmer ^ kmer >> 28;
	kmer = (kmer + (kmer << 31)) & mask;
	return kmer;
}

//This function calculates the numerical representation of a k-mer
uint64_t calcKmerNb(const string& kmer);

//This function builds the sketch of a sequence using a given hash threshold
PairSketch buildSketch(const string& s, const uint64_t& ht);

//This function builds a FracMinHash sketch of a sequence using a given hash threshold. K-mers occurring on the given black list are
//ignored
const Sketch buildSketch(const string& seq, const uint32_t& k, const double& hFrac, const unordered_map<uint64_t, char>& bLstmers);

//This function builds a minimap2 sketch of a sequence by querying from a prebuilt minimap index
const Sketch buildMiniSketch(const string& seq, const mm_idx_t *pidx);

//This function dereplicates elements from sketches sharing the same hash value (only the first one is kept). WARNING: It does not preserve elements'
//initial ordering
void remDuplHshs(PairSketch& sk);

#endif