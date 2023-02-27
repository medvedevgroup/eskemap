#include "Sketch.h"
#include <math.h>

//This function calculates the numerical representation of a k-mer
uint64_t calcKmerNb(const string& kmer){
	uint64_t kmerNb = 0;

	for(string::const_iterator c = kmer.begin(); c != kmer.end(); ++c) kmerNb = (kmerNb << 2) + getBaseNb(*c);

	return kmerNb;
}

//This function builds the sketch of a sequence
PairSketch buildSketch(const string& s){
	const uint64_t mask = pow(4, K) - 1;
	uint64_t kmerHash;
	PairSketch sk;

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

//This function builds a FracMinHash sketch of a sequence using a given hash threshold. K-mers occurring on the given black list are
//ignored
const Sketch buildSketch(const string& seq, const uint32_t& k, const double& hFrac, const unordered_map<uint64_t, char>& bLstmers){
	const uint64_t mask = pow(ALPHABET_SIZE, k) - 1;
	//Calculate maximum hash value in sketch
	const uint64_t maxHash = hFrac * mask;
	uint64_t kmerHash;
	Sketch sk;

	//If the sequence is smaller than k we are done
	if(seq.length() < k){
		cerr << "WARNING: Length of input sequence " << seq << " is smaller than k (k=" << k << ")" << endl;

		return sk;
	}

	//Reserve as much space as is approximately needed to store the sketch (which hopefully saves some time)
	sk.reserve(seq.length() * hFrac);

	//Iterate over k-mer starting positions in sequence
	for(uint32_t i = 0; i < seq.length() - k + 1; ++i){
		//Calculate numerical k-mer representation and its hash
		kmerHash = getHash(calcKmerNb(seq.substr(i, k)), mask);

		//Check if hash value is small enough to be kept
		if(kmerHash <= maxHash && !bLstmers.contains(kmerHash)) sk.push_back(kmerHash);
	}

	//Resize sketch (just for case we have allocated way too much memory)
	sk.shrink_to_fit();

	return sk;
}

//This function builds a minimap2 sketch of a sequence by querying from a prebuilt minimap index; this function is influenced by the
// code of "The minimizer Jaccard estimator is biased and inconsistent." from Belbasi et al. (function 
//"winnowed_minimizers_linear(perm,windowSize)" in file "winnowed_minimizers.py").
const Sketch buildMiniSketch(const string& seq, const uint32_t& k, const uint32_t& w, const unordered_map<uint64_t, char>& blmers){
	uint32_t lastIdx = UINT32_MAX;
	int32_t windowBorder;
	const uint64_t mask = pow(ALPHABET_SIZE, k) - 1;
	uint64_t kmerHash, kmerBitSeq, revKmerBitSeq;
	//This stores pairs of k-mer starting positions and their hashes within the current window
	vector<pair<int32_t, uint64_t>> windowKmers;
	Sketch sk;

	//If the sequence is smaller than k we are done
	if(seq.length() < k){
		cerr << "WARNING: Length of input sequence " << seq << " is smaller than k (k=" << k << ")" << endl;

		return sk;
	}

	//Reserve as much space as is approximately needed to store the sketch (which hopefully saves some time)
	sk.reserve(seq.length() * 0.03);

	//Iterate over k-mer starting positions in sequence
	for(int32_t i = 0; i < seq.length() - k + 1; ++i){
		kmerBitSeq = calcKmerNb(seq.substr(i, k));
		revKmerBitSeq = calcKmerNb(revComp(seq.substr(i, k)));
		windowBorder = i - (w - 1);

		//We do not consider k-mers which are their own reverse complement
		if(kmerBitSeq == revKmerBitSeq) continue;

		//Depending on which is lexicographically smaller, we consider either the k-mer or its reverse complement
		if(kmerBitSeq < revKmerBitSeq){
			//Calculate k-mer's hash
			kmerHash = getHash(kmerBitSeq, mask);
		} else{
			//Calculate k-mer's hash
			kmerHash = getHash(revKmerBitSeq, mask);
		}
		
		//Remove all pairs with a larger hash from the back
		while(!windowKmers.empty() && windowKmers.back().second > kmerHash) windowKmers.pop_back();

		//Add current k-mer hash to windowKmers
		windowKmers.push_back(make_pair(i, kmerHash));

		//Remove pairs of k-mers which are no longer inside the window
		while(!windowKmers.empty() && windowKmers.front().first < windowBorder) windowKmers.erase(windowKmers.begin());

		//Choose a minimizer as soon as we have the first full window of k-mers and make sure we do the same minimizer a second time
		if(windowBorder >= 0 && !windowKmers.empty()){
			if(lastIdx != windowKmers.front().first && !blmers.contains(windowKmers.front().second)){
				lastIdx = windowKmers.front().first;
				sk.push_back(windowKmers.front().second);
			}

			//If the same k-mer appears several times inside the window and it has the smallest hash we want to save all occurrences
			while(windowKmers.size() > 1 && windowKmers.front().second == windowKmers[1].second){
				windowKmers.erase(windowKmers.begin());
				lastIdx = windowKmers.front().first;

				//Blacklisted k-mers are not added
				if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front().second);
			}
		}
	}

	//In case we have never seen a full window of k-mers take the one with the smallest hash seen for the sketch
	if(windowBorder < 0 && !windowKmers.empty()){
		//Blacklisted k-mers are not added
		if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front().second);

		//If the same k-mer appears several times inside the window and it has the smallest hash we want to save all occurrences
		while(windowKmers.size() > 1 && windowKmers.front().second == windowKmers[1].second){
			windowKmers.erase(windowKmers.begin());

			//Blacklisted k-mers are not added
			if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front().second);
		}
	}

	//Resize sketch (just for case we have allocated way too much memory)
	sk.shrink_to_fit();

	return sk;
}

//This function dereplicates elements from sketches sharing the same hash value (only the first one is kept). WARNING: It does not preserve elements'
//initial ordering
void remDuplHshs(PairSketch& sk){
	PairSketch::iterator i, j;

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

//This function calculates the reverse complement of a DNA sequence
string revComp(const string &seq){
	//The result string
	string revSeq;

	//Go through the query from the end to the beginning
	for(int32_t i = seq.length() - 1; i >= 0; --i){
		//Check which base we are dealing with and append its complement
		switch(seq[i]){
			case NUCL_BASE_A:
				revSeq += CMPL_BASE_A;
				break;
			case NUCL_BASE_C:
				revSeq += CMPL_BASE_C;
				break;
			case NUCL_BASE_G:
				revSeq += CMPL_BASE_G;
				break;
			case NUCL_BASE_T:
				revSeq += CMPL_BASE_T;
				break;
			default:
				revSeq += "N";
				break;
		}
	}

	return revSeq;
}
