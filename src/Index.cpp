#include "Index.h"
#include "algorithm"

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
