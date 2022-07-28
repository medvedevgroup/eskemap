#include <unordered_map>

#include "Thomology.h"

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const Sketch& skT, const uint32_t& cw, const uint32_t& uw, const int32_t& t, 
	const bool& norm){
	uint32_t k;
	//This variable stores delta_occ
	int32_t delta
	//An array to store a mapping between columns and positions in the text they belong to
	vector<uint32_t> col2pos;
	vector<Thomology> res;
	//The score matrix
	vector<vector<int32_t>> scores;
	//Initialize abcs (we assume there are no duplicates in p)
	unordered_map<uint64_t, uint32_t> abcs(skP.size());
	//Initialize occt (again, we assume there are no duplicates in p)
	unordered_map<uint64_t, vector<uint32_t>> occt(skP.size());

	//Fill abc
	for(Sketch::const_iterator i = skP.begin(); i != skP.end(); ++i){
		if(abcs.contains(*i)){
			++abcs[*i];
		} else{
			abcs[i] = 1;
		}
	}

	//Fill score matrix
	for(uint32_t j = 0; j < skT.size(); ++j){
		//If hash does not exist in pattern we can skip this position
		if(!abcs.contains(skT[j])) continue;

		//Save this columns position in the text
		col2pos.push_back(j);
		//We don't want to search k_min from the beginning again each time
		k = 0;

		//Iterate over start positions
		for(vector<uint32_t>::const_iterator i = col2pos.begin(); i != col2pos.end(); ++i){
			//Check if current hash occurred already before
			if(occt.contains(skT[j])){
				//Walk through vector at current position
				while(k < occt[skT[j]].size()){
					//Check if we have found k_min already
					if(*i < occt[skT][k]){
						//Calculate delta_occ
						delta = abcs[skT[j]] - occt[skT[j]].size() + k;
						break;//Maybe better use goto here?
					}

					++k;
				}


			} else{
				occt[skT] = {j};
				//Calculate delta_occ
				delta = abcs[skT[j]] - 1;//At this point we know in principle that delta_occ is non-negative
			}
			
		}
	}

	//Find maximum t-homologies

	return res;
}