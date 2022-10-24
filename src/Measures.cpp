#include <vector>
#include <algorithm>

#include "Measures.h"
#include "Sketch.h"

//This function calculates the "set intersection" similarity score
const int32_t calcIntersecScore(PairSketch& skA, PairSketch& skB){
	bool found;
	uint32_t nbShared = 0, nbUniqA = 0;

	//Remove duplicated hashes from sketches
	remDuplHshs(skA);
	remDuplHshs(skB);

	//Iterate over first sketch
	for(PairSketch::const_iterator i = skA.begin(); i != skA.end(); ++i){
		//Reset flag
		found = false;

		//Iterate over second sketch
		for(PairSketch::const_iterator j = skB.begin(); j != skB.end(); ++j){
			//Check if both current elements share the same hash value
			if(i->second == j->second){
				++nbShared;
				found = true;
				break;
			}
		}

		//Increment counter if current hash value could not be found in other sketch
		if(!found) ++nbUniqA;
	}

	//shared - uniqA - uniqB = 2 * shared - uniqA - size(B)
	return 2 * nbShared - nbUniqA - skB.size();
}

//This function calculates the "maximum aligned hashes" similarity score
const int32_t calcAlgnHshsScore(const PairSketch& skA, const PairSketch& skB, const bool& consOffs){
	//The matrix rows
	vector<int32_t> lst, cur;
	//Sketch iterators
	PairSketch::const_iterator ia, ib = skB.begin();

	//Initialize first row
	for(int32_t j = 0; j <= skA.size(); ++j) lst.push_back(-j);

	//Calculate all remaining values
	for(int32_t i = 1; i <= skB.size(); ++i){
		//Initialize iterator of sketch A
		ia = skA.begin();

		//Calculate current row
		for(int32_t j = 0; j <= skA.size(); ++j){
			//Check if the current element is the first in a row
			if(j == 0){
				//For the first element we do not need to compare anything
				cur.push_back(-i);
			} else{
				//Calculate the current element
				cur.push_back(max(lst[j-1] + (ia->second == ib->second ? 1 : -1), max(lst[j] - 1, cur[j-1] - 1)));
				//Increment iterator
				++ia;
			}
		}

		//Make current row the last one
		lst = cur;
		//Reset current row
		cur.clear();
		//Increment iterator
		++ib;
	}

	return lst.back();
}
