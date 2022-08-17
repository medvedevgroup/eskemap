#include <unordered_map>

#include "Thomology.h"

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const Sketch& skT, const uint32_t& cw, const uint32_t& uw, const int32_t& t){
	//Some counter variables
	uint32_t i, j, k, occ;
	//The maximum threshold to compare against
	int32_t maxThres;
	//An array to store a mapping between columns and positions in the text they belong to
	vector<uint32_t> col2pos;
	vector<Thomology> res;
	//The score matrix
	vector<vector<int32_t>> scores;
	//A sketch to store a reduced representation of the text
	Sketch L;
	//An iterator to iterate over L
	Sketch::const_iterator fSkIt;
	//An iterator to reversely iterate over L
	Sketch::const_reverse_iterator rSkIt;
	//An iterator to iterate over lists in pos
	vector<uint32_t>::const_iterator posIt;
	//A list to store scores of maximum t-homologies along with their start position
	list<pair<uint32_t, int32_t>> maxScores;
	//An iterator for the maximum scores list
	list<pair<uint32_t, int32_t>>::iterator li;
	//Initialize occp (we assume there are no duplicates in p)
	unordered_map<uint64_t, uint32_t> occp(skP.size());
	//Initialize pos (again, we assume there are no duplicates in t)
	unordered_map<uint64_t, vector<uint32_t>> pos(skP.size());

	//Fill occp
	for(fSkIt = skP.begin(); fSkIt != skP.end(); ++fSkIt){
		if(occp.contains(*fSkIt)){
			++occp[*fSkIt];
		} else{
			occp[*fSkIt] = 1;
		}
	}

	//Set position counter
	j = 0;
	//Reserve some space for L and col2pos
	L.reserve(P_MULTIPLICITY * skP.size());
	col2pos.reserve(P_MULTIPLICITY * skP.size());

	//Reduce the text on matching positions
	for(fSkIt = skT.begin(); fSkIt != skT.end(); ++fSkIt){
		//Check if k-mer occurs in pattern
		if(occp.contains(*fSkIt)){
			//Store k-mer in L
			L.push_back(*fSkIt);
			//Store k-mer's position in t
			col2pos.push_back(j);

			//Testing
			if(*fSkIt == 74701964) cout << "findThoms: j: " << j << endl;
		}

		++j;
	}

	//From now on L and col2pos will keep their sizes
	L.shrink_to_fit();
	col2pos.shrink_to_fit();
	//Reset position counter
	j = 0;

	//Fill score matrix
	for(fSkIt = L.begin(); fSkIt != L.end(); ++fSkIt){
		//Add new column in score matrix...
		scores.push_back(vector<int32_t>());
		//...and reserve enough space
		scores.back().reserve(j + 1);

		//Check if current hash occurred already before
		if(pos.contains(*fSkIt)){
			//Add current position as occurrence for this hash
			pos[*fSkIt].push_back(j);
		} else{
			pos[*fSkIt] = {j};
		}

		//We don't want to search k_min from the beginning again each time
		k = 0;
		posIt = pos[*fSkIt].begin();

		//Iterate over start positions
		for(i = 0; i <= j; ++i){
			//Walk through vector at current position
			while(posIt != pos[*fSkIt].end()){
				//Check if we have found k_min already
				if(i <= *posIt){
					//Calculate occ(t[j], t[i, j-1])
					occ = pos[*fSkIt].size() - k - 1;
					break;
				}

				++k;
				++posIt;
			}

			//Check if we are dealing with the base case
			if(i == j){
				scores.back().push_back(cw - uw * (skP.size() - 1));
			} else if(occ < occp[*fSkIt]){//Discriminate between cases
				scores.back().push_back(scores[j - 1][i] + cw + uw - uw * (col2pos[j] - col2pos[j - 1] - 1));
			} else{
				scores.back().push_back(scores[j - 1][i] - uw * (col2pos[j] - col2pos[j - 1]));
			}
		}

		++j;
	}

	//Get a reverse iterator to iterate over L
	rSkIt = L.rbegin();
	//Get a counter for the column that we are at
	j = scores.size() - 1;

	//Find maximum t-homologies
	for(vector<vector<int32_t>>::const_reverse_iterator colRit = scores.rbegin(); colRit != scores.rend(); ++colRit){
		//Reset row counter
		i = 0;
		//Start in the beginning of the list
		li = maxScores.begin();
		//Set the threshold that a maximum t-homology has to exceed
		maxThres = t - 1;
		//Get a forward iterator to iterate over L
		fSkIt = L.begin();

		//Walk through column from top to bottom
		for(vector<int32_t>::const_iterator rowIt = colRit->begin(); rowIt != colRit->end(); ++rowIt){
			//If the substring is a t-homology and its first and last hash is identical this hash needs to occur at least twice inside the pattern
			if(i != j && *fSkIt == *rSkIt && occp[*fSkIt] < 2){
				++i;
				++fSkIt;
				continue;
			}

			//Walk through the list to consider all relevant maximums seen so far
			while(li != maxScores.end()){
				//Maximums found in rows > i are irrelevant at this point
				if(i < li->first) break;

				//Update maximum to compare with
				maxThres = max(maxThres, li->second);
				//Walk on
				++li;
			}

			//Check if score is high enough
			if(*rowIt > maxThres){
				//Testing
				cout << "findThoms: *fSkIt: " << *fSkIt << " *rSkIt: " << *rSkIt << endl;

				//Add t-homology to results
				res.push_back(make_tuple(col2pos[i], col2pos[j], *rowIt));
				//Add score to list with maximum scores
				maxScores.insert(li, make_pair(i, *rowIt));
				//Update maximum to comare with
				maxThres = *rowIt;
			}

			++i;
			++fSkIt;
		}

		++rSkIt;
		--j;
	}

	return res;
}
