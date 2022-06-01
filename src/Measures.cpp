#include "Measures.h"
#include "Sketch.h"

//This function calculates the "set intersection" similarity score
const int32_t calcIntersecScore(Sketch& skA, Sketch& skB){
	bool found;
	uint32_t nbShared = 0, nbUniqA = 0;

	//Remove duplicated hashes from sketches
	remDuplHshs(skA);
	remDuplHshs(skB);

	//Iterate over first sketch
	for(Sketch::const_iterator i = skA.begin(); i != skA.end(); ++i){
		//Reset flag
		found = false;

		//Iterate over second sketch
		for(Sketch::const_iterator j = skB.begin(); j != skB.end(); ++j){
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
