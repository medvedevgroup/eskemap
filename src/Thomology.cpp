#include <unordered_map>

#include "Thomology.h"

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const Sketch& skT, const uint32_t& cw, const uint32_t& uw, const int32_t& t){
	uint32_t i, j = 0, k, occ;
	//An array to store a mapping between columns and positions in the text they belong to
	vector<uint32_t> col2pos(P_MULTIPLICITY * skP.size());
	vector<Thomology> res;
	//The score matrix
	vector<vector<int32_t>> scores;
	//A sketch to store a reduced representation of the text
	Sketch L(P_MULTIPLICITY * skP.size());
	//An iterator to iterate over lists in pos
	vector<uint32_t>::const_iterator posIt;
	//An array storing relevant maximum scores for each row
	unordered_map<uint32_t, int32_t> maxScores;
	//Initialize occp (we assume there are no duplicates in p)
	unordered_map<uint64_t, uint32_t> occp(skP.size());
	//Initialize pos (again, we assume there are no duplicates in t)
	unordered_map<uint64_t, vector<uint32_t>> pos(skP.size());

	//Fill occp
	for(Sketch::const_iterator i = skP.begin(); i != skP.end(); ++i){
		if(occp.contains(*i)){
			++occp[*i];
		} else{
			occp[*i] = 1;
		}
	}

	//Reduce the text on matching positions
	for(Sketch::const_iterator i = skT.begin(); i != skT.end(); ++i){
		//Check if k-mer occurs in pattern
		if(occp.contains(*i)){
			//Store k-mer in L
			L.push_back(*i);
			//Store k-mer's position in t
			col2pos.push_back(j);
		}
		++j;
	}

	//From now on L and col2pos will keep their sizes
	L.shrink_to_fit();
	col2pos.shrink_to_fit();
	//Reset position counter
	j = 0;

	//Fill score matrix
	for(Sketch::const_iterator colIt = L.begin(); colIt != L.end(); ++colIt){
		//Add new column in score matrix...
		scores.push_back(vector<int32_t>());
		//...and reserve enough space
		scores.back().reserve(j + 1);
		//Counter to know which line we are currently at
		i = 0;

		//Check if current hash occurred already before
		if(pos.contains(*colIt)){
			//Add current position as occurrence for this hash
			pos[*colIt].push_back(j);
		} else{
			pos[*colIt] = {j};
		}

		//We don't want to search k_min from the beginning again each time
		k = 0;
		posIt = pos[*colIt].begin();

		//Iterate over start positions
		for(Sketch::const_iterator rowIt = L.begin(); i <= j; ++rowIt){
			//Walk through vector at current position
			while(posIt != pos[*colIt].end()){
				//Check if we have found k_min already
				if(i <= *posIt){
					//Calculate occ(t[j], t[i, j-1])
					occ = pos[*colIt].size() - k - 1;
					break;
				}

				++k;
				++posIt;
			}

			//Check if we are dealing with the base case
			if(i == j){
				scores.back().push_back(cw - uw * (skP.size() - 1));
			} else if(occ < occp[*colIt]){//Discriminate between cases
				scores.back().push_back(scores[j - 1][i] + cw + uw - uw * (col2pos[j] - col2pos[j - 1] - 1));
			} else{
				scores.back().push_back(scores[j - 1][i] - uw * (col2pos[j] - col2pos[j - 1]));
			}

			++i;
		}

		++j;
	}

	//Find maximum t-homologies
	for(vector<vector<int32_t>>::const_reverse_iterator colRit = scores.rbegin(); colRit != scores.end(); ++colRit){
		//Walk through column from top to bottom
		for(vector<int32_t>::const_iterator rit = colRit->begin(); colRit->end(); ++rit){
			//Check if we had a maximum for previous row
			if(maxScores.contains(i - 1)){
				//Check if we had a maximum for this row already, too
				if(maxScores.contains(i)){
					maxScores[i] = max(maxScores[i], maxScores[i - 1]);
				} else{
					maxScores[i] = maxScores[i - 1];
				}
			}

			//If the substing is a t-homology and its first and last hash is identical this hash needs to occur at least twice inside the pattern
			if(skT[col2pos[i]] == skT[col2pos[j]] && occp[skT[col2pos[j]]] < 2) continue;

			//Check if we have already seen a relevant maximum
			if(maxScores.contains(i)){
				//Check if current substring is maximal
				if(scores[j][i] > maxScores[i]){
					//Update maximum
					maxScores[i] = scores[j][i];
					//Add t-homology to results
					res.push_back(make_tuple(col2pos[i], col2pos[j], scores[j][i]));
				} else if(scores[j][i] >= t){
					//Save new maximum
					maxScores[i] = scores[j][i];
					//Add t-homology to results
					res.push_back(make_tuple(col2pos[i], col2pos[j], scores[j][i]));
				}
			}
		}
	}

	return res;
}