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

//This function builds a minimap2 sketch of a sequence by querying from a prebuilt minimap index
const Sketch buildMiniSketch(const string& seq, const string& id, const mm_idx_t *pidx){
	int nHits;
	uint32_t sIdMini;
	const uint64_t mask = pow(ALPHABET_SIZE, pidx->k) - 1;
	uint64_t kmerHash;
	const mm_idx_seq_t *seqFeats = pidx->seq;
	Sketch sk;

	//If the sequence is smaller than k we are done
	if(seq.length() < pidx->k){
		cerr << "WARNING: Length of input sequence " << seq << " is smaller than k (k=" << pidx->k << ")" << endl;

		return sk;
	}

	//Reserve as much space as is approximately needed to store the sketch (which hopefully saves some time)
	sk.reserve(seq.length() * 0.03);

	//Find index internal reference id of our sequence
	for(uint32_t i = 0; i < pidx->n_seq; ++i){
		if(!strcmp((seqFeats+i)->name, id.c_str())){
			sIdMini = i;

			//Testing
			// cout << "buildMiniSketch: read " << id << " has the internal id " << sIdMini << endl;

			break;
		}
	}

	//Iterate over k-mer starting positions in sequence
	for(uint32_t i = 0; i < seq.length() - pidx->k + 1; ++i){
		//Calculate numerical k-mer representation and its hash
		kmerHash = getHash(calcKmerNb(seq.substr(i, pidx->k)), mask);
		//Query current hash in index
		const uint64_t *idx_p = mm_idx_get(pidx, kmerHash, &nHits);

		//Check if hash could be found
		if(nHits > 0){
			//Iterate over all occurrences
		    for(uint32_t j = 0; j < nHits; ++j){
		    	if((uint32_t) ((*idx_p)>>32) == sIdMini){
		    		sk.push_back(kmerHash);
		    		break;
		    	}

		        //Move to next occurrence
		        idx_p++;
	        }
	    }

	    //Testing
	    if(kmerHash == 26942362073) cout << "Interesting hash queried in forward direction" << endl;
	    uint64_t lastKmerHash = kmerHash;

	    //Do the same for the reverse complement
	    kmerHash = getHash(calcKmerNb(revComp(seq.substr(i, pidx->k))), mask);

	    //Testing
	    if(kmerHash == 26942362073) cout << "Interesting hash queried in backward direction" << endl;
	    if(kmerHash == lastKmerHash) cerr << "buildMiniSketch: We hash for k-mer and its reverse complement are the same" << endl;

		const uint64_t *idx_q = mm_idx_get(pidx, kmerHash, &nHits);

		//Testing
	    // if(kmerHash == 16105810989 && nHits > 0) cout << "buildMiniSketch: ...and could be found" << endl;

		//Check if hash could be found
		if(nHits > 0){
			//Iterate over all occurrences
		    for(uint32_t j = 0; j < nHits; ++j){
		    	if((uint32_t) ((*idx_q)>>32) == sIdMini){
		    		sk.push_back(kmerHash);
		    		break;
		    	}

		        //Move to next occurrence
		        idx_q++;
	        }
	    }

		//Testing
		// cout << "buildMiniSketch: Read id is " << ((*idx_p)>>32) << endl;
		// cout << "buildMiniSketch: Searching for k-mer hash " << kmerHash << endl;
		// cout << "buildMiniSketch: i: " << i << 
	}

	//Testing
	// const mm_idx_seq_t *s = pidx->seq;
	// cout << s->name << endl;
	// cout << (s+1)->name << endl;

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
