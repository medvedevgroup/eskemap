#ifndef INDEX_HPP
#define INDEX_HPP

#include <utility>

//A compare function to sort elements in a hash position vector
inline const bool smHshnPos(const pair<uint64_t, uint32_t>&hpa, const pair<uint64_t, uint32_t>& hpb){//TODO: This function still needs to be tested!
	//Testing
	if(hpa.first == hpb.first){
		cout << "1 Option " << (hpa.second < hpb.second ? "1" : "2") << endl;
	} else{
		cout << "2 Option " << (hpa.second < hpb.second ? "1" : "2") << endl;
	}

	return hpa.first != hpb.first ? hpa.first < hpb.first : hpa.second < hpb.second;
}

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