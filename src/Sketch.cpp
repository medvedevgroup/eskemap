#include "Sketch.h"
#include <math.h>

//This function calculates the numerical representation of a k-mer
uint64_t calcKmerNb(const string& kmer){
	uint64_t kmerNb = 0;

	for(string::const_iterator c = kmer.begin(); c != kmer.end(); ++c) kmerNb = (kmerNb << 2) + getBaseNb(*c);

	return kmerNb;
}

//This function builds the sketch of a sequence
Sketch buildSketch(const string& s){
	const uint64_t mask = pow(4, K) - 1;
	uint64_t kmerHash;
	list<pair<uint32_t, uint64_t>> sk;

	//If the sequence is smaller than k we are done
	if(s.length() < K){
		cerr << "WARNING: Length of input sequence " << s << " is smaller than k (k=" << K << ")" << endl;

		return sk;
	}

	//Iterate over k-mer starting positions in sequence
	for(uint32_t i = 0; i < s.length() - K + 1; ++i){
		//Calculate numerical k-mer representation and its hash
		kmerHash = getHash(calcKmerNb(s.substr(i, K)), mask);

		//Check if hash value is small enough to be kept
		if(kmerHash <= MAX_HASH) sk.push_back(make_pair(i, kmerHash));
	}

	return sk;
}

//This function dereplicates elements from sketches sharing the same hash value (only the first one is kept). WARNING: It does not preserve elements'
//initial ordering
void remDuplHshs(Sketch& sk){
	Sketch::iterator i, j;

	//If the sketch consists of only one element we are done
	if(sk.size() < 2) return;

	//Sort elements according to their hash values
	sk.sort(smHshsFrst);
	//Iterate over elements
	i = sk.begin();
	j = sk.begin();
	//j is one element in front of i
	++j;

	while(j != sk.end()){
		//If elements share the same hash value delete the second
		if(i->second == j->second){
			j = sk.erase(j);
		} else{
			//Move to next elements
			++i;
			++j;
		}
	}
}
