#include "Thomology.h"
#include "Index.h"

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const mm_idx_t *tidx, const uint32_t& cw, 
	const uint32_t& uw, const int32_t& t){
	//Some counter variables
	uint32_t i, j, k, occ;
	//The maximum threshold to compare against
	int32_t maxThres;
	vector<Thomology> res;
	//The score matrix
	vector<vector<int32_t>> scores;
	//A sketch interator
	Sketch::const_iterator fSkIt;
	//Hash position array of all hashes and their positions in the text sketch which also appear inside the pattern
	vector<pair<uint64_t, uint32_t>> L;
	//An iterator to iterate over L
	vector<pair<uint64_t, uint32_t>>::const_iterator fLit;
	//An iterator to reversely iterate over L
	vector<pair<uint64_t, uint32_t>>::const_reverse_iterator rLit;
	//An iterator to iterate over lists in pos
	vector<uint32_t>::const_iterator posIt;
	//A list to store scores of maximum t-homologies along with their start position
	list<pair<uint32_t, int32_t>> maxScores;
	//An iterator for the maximum scores list
	list<pair<uint32_t, int32_t>>::iterator li;
	//Initialize occp
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

	//Testing
	// cout << "findThoms: Filling of occp done" << endl;

	//Generate L
	L = genL(occp, tidx);

	//Testing
	// cout << "findThoms: Generated L" << endl;
	// cout << "findThoms: Its size is " << L.size() << endl;

	//Set position counter
	j = 0;

	//Fill score matrix
	for(fLit = L.begin(); fLit != L.end(); ++fLit){
		//Add new column in score matrix...
		scores.push_back(vector<int32_t>());
		//...and reserve enough space
		scores.back().reserve(j + 1);

		//Check if current hash occurred already before//TODO: Can we speed this up by using the index?
		if(pos.contains(fLit->first)){
			//Add current position as occurrence for this hash
			pos[fLit->first].push_back(j);
		} else{
			pos[fLit->first] = {j};
		}

		//We don't want to search k_min from the beginning again each time
		k = 0;
		posIt = pos[fLit->first].begin();

		//Iterate over start positions
		for(i = 0; i <= j; ++i){
			//Walk through vector at current position
			while(posIt != pos[fLit->first].end()){
				//Check if we have found k_min already
				if(i <= *posIt){
					//Calculate occ(t[j], t[i, j-1])
					occ = pos[fLit->first].size() - k - 1;
					break;
				}

				++k;
				++posIt;
			}

			//Check if we are dealing with the base case
			if(i == j){
				scores.back().push_back(cw - uw * (skP.size() - 1));
			} else if(occ < occp[fLit->first]){//Discriminate between cases
				scores.back().push_back(scores[j - 1][i] + cw + uw - uw * (fLit->second - L[j - 1].second - 1));
			} else{
				scores.back().push_back(scores[j - 1][i] - uw * (fLit->second - L[j - 1].second));
			}
		}

		++j;
	}

	//Testing
	// cout << "findThoms: Scores calculated" << endl;
	// cout << "findThoms: Size of score matrix: " << scores.size() << "^2" << endl;

	//Get a reverse iterator to iterate over L
	rLit = L.rbegin();
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
		fLit = L.begin();

		//Walk through column from top to bottom
		for(vector<int32_t>::const_iterator rowIt = colRit->begin(); rowIt != colRit->end(); ++rowIt){
			//If the substring is a t-homology and its first and last hashes are identical this hash needs to occur at least twice 
			//inside the pattern
			if(i != j && fLit->first == rLit->first && occp[fLit->first] < 2){
				++i;
				++fLit;
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
				//Add t-homology to results
				res.push_back(make_tuple(L[i].second, L[j].second, *rowIt));
				//Add score to list with maximum scores
				maxScores.insert(li, make_pair(i, *rowIt));
				//Update maximum to comare with
				maxThres = *rowIt;
			}

			++i;
			++fLit;
		}

		++rLit;
		--j;
	}

	return res;
}
