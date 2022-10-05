#include "Index.h"
#include "../../software/minimap2/index.c"

//Index constructor
Index::Index(const Sketch& sk){
	uint32_t i = 0, currPos = 0, start, end;
	uint64_t lastHash;
	vector<pair<uint64_t, uint32_t>>::const_iterator e;

	//Reserve enough space for vector
	this->pos.reserve(sk.size());

	//Store hash position pairs in pos
	for(Sketch::const_iterator h = sk.begin(); h != sk.end(); ++h, ++i) this->pos.push_back(make_pair(*h, i));

	//Resize pos after filling
	this->pos.shrink_to_fit();
	//Sort elements
	sort(this->pos.begin(), this->pos.end(), smHshnPos);
	//Fill hash table
	e = this->pos.begin();
	lastHash = e->first;
	start = 0;
	end = 0;

	while(e != this->pos.end()){
		if(lastHash != e->first){
			this->pints[lastHash] = make_pair(start, end);
			start = currPos;
			end = currPos;
			lastHash = e->first;
		}

		++currPos;
		++end;
		++e;
	}

	this->pints[lastHash] = make_pair(start, end);
}

//This function generates the L array need for alpha-homology detection
const vector<pair<uint64_t, uint32_t>> genL(const Sketch& pSk, const mm_idx_t *tidx, const uint32_t& k){
	uint32_t i;
	//Number of times a hash can be found inside the text sketch
	int nHits;
	//An iterator to iterate over the sketch
	Sketch::const_iterator sI;
	//The L array
	vector<pair<uint64_t, uint32_t>> L;

	//Reserve some space for L
	L.reserve(P_MULTIPLICITY * pSk.size());

	//Query index for all hashes occurring in the read sketch
	for(sI = pSk.begin(); sI != pSk.end(); ++sI){
		//Query current hash in index
		const uint64_t *idx_p = mm_idx_get(tidx, *sI, &nHits);

		//Check if hash could be found
		if(nHits > 0){
		    //Iterate over all occurrences
		    for(i = 0; i < nHits; ++i){
		    	L.push_back(make_pair(*sI, (((uint32_t)(*idx_p))>>1) - k + 1));
		        //Move to next occurrence
		        idx_p++;
	        }
		}
	}

	//Sort L by ascending positions in text
	sort(L.begin(), L.end(), smHshnPos);

	return L;
}
