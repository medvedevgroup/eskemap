#include "Sketch.h"
#include <math.h>

//This function calculates the numerical representation of a k-mer
uint64_t calcKmerNb(const std::string& kmer){
	uint64_t kmerNb = 0;

	for(std::string::const_iterator c = kmer.begin(); c != kmer.end(); ++c) kmerNb = (kmerNb << 2) + getBaseNb(*c);

	return kmerNb;
}

//This function builds the sketch of a sequence
Sketch buildSketch(const std::string& s){
	const uint64_t mask = pow(4, K) - 1;
	uint64_t kmerHash;
	std::list<std::pair<uint32_t, uint64_t>> sk;

	//Iterate over k-mer starting positions in sequence
	for(uint32_t i = 0; i < s.length() - K + 1; ++i){
		//Calculate numerical k-mer representation and its hash
		kmerHash = getHash(calcKmerNb(s.substr(i, K)), mask);//This function still needs to be implemented!

		//Check if hash value is small enough to be kept
		if(kmerHash <= MAX_HASH) sk.push_back(std::make_pair(i, kmerHash));
	}

	return sk;
}