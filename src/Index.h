#ifndef INDEX_HPP
#define INDEX_HPP

#include "Sketch.h"
#include "../../software/minimap2/minimap.h"

#define INDEX_DEFAULT_DUMP_FILE "indexDump.idx"

//A compare function to sort elements in a hash position vector primarily by hash and secondarily by position in ascending order
inline const bool smHshnPos(const pair<uint64_t, uint32_t>& hpa, const pair<uint64_t, uint32_t>& hpb){
	return hpa.first != hpb.first ? hpa.first < hpb.first : hpa.second < hpb.second;
}

//A compare function to sort elements in a hash position vector by ascending positions only
inline const bool smPos(const pair<uint64_t, uint32_t>& hpa, const pair<uint64_t, uint32_t>& hpb){
	return hpa.second < hpb.second;
}

//This function generates the L array need for alpha-homology detection
const vector<pair<uint64_t, uint32_t>> genL(const unordered_map<uint64_t, uint32_t>& phshs, const mm_idx_t *tidx);

class Index {

public:

		Index(const Sketch& sk);

protected:

		//A hash table storing occurrence intervals of k-mers in pos (hash -> [start, end))
		unordered_map<uint64_t, pair<uint32_t, uint32_t>> pints;
		//Vector of k-mer hash position pairs
		vector<pair<uint64_t, uint32_t>> pos;
};

#endif