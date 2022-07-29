#include <unordered_map>

#include "Thomology.h"

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const Sketch& skT, const uint32_t& cw, const uint32_t& uw, const int32_t& t){
	uint32_t k;
	//This variable stores delta_occ
	int32_t delta;
	//An array to store a mapping between columns and positions in the text they belong to
	vector<uint32_t> col2pos;
	vector<Thomology> res;
	//The score matrix
	vector<vector<int32_t>> scores;
	//An array storing relevant maximum scores for each row
	unordered_map<uint32_t, int32_t> maxScores;
	//Initialize abcs (we assume there are no duplicates in p)
	unordered_map<uint64_t, uint32_t> abcs(skP.size());
	//Initialize occt (again, we assume there are no duplicates in p)
	unordered_map<uint64_t, vector<uint32_t>> occt(skP.size());

	//Fill abc
	for(Sketch::const_iterator i = skP.begin(); i != skP.end(); ++i){
		if(abcs.contains(*i)){
			++abcs[*i];
		} else{
			abcs[*i] = 1;
		}
	}

	//Fill score matrix
	for(uint32_t j = 0; j < skT.size(); ++j){
		//If hash does not exist in pattern we can skip this position
		if(!abcs.contains(skT[j])) continue;

		//Add new column in score matrix...
		scores.push_back(vector<int32_t>());
		//...and reserve enough space
		scores.back().reserve(j + 1);
		//Save this columns position in the text
		col2pos.push_back(j);
		//We don't want to search k_min from the beginning again each time
		k = 0;

		//Check if current hash occurred already before
		if(occt.contains(skT[j])){
			//Add current position as occurrence for this hash
			occt[skT[j]].push_back(j);
		} else{
			occt[skT[j]] = {j};
		}

		//Iterate over start positions (i corresponds to the current row in the score matrix)
		for(uint32_t i = 0; i < col2pos.size(); ++i){
			//Walk through vector at current position
			while(k < occt[skT[j]].size()){
				//Check if we have found k_min already
				if(col2pos[i] <= occt[skT[j]][k]){
					//Calculate delta_occ
					delta = abcs[skT[j]] - occt[skT[j]].size() + k;
					break;
				}

				++k;
			}

			//Check if substring of score to be calculated starts at position j
			if(col2pos[i] == j){
				//At this point we know in principle that delta_occ is non-negative
				scores.back().push_back(2 * cw - (skP.size() - 1) * uw);
			} else if(delta < 0){//Check if delta increased
				scores.back().push_back(scores[j - 1][i] - uw);
			} else{
				scores.back().push_back(scores[j - 1][i] + 2 * cw + uw);
			}
		}
	}

	//Find maximum t-homologies (This time j is position in score matrix not in text)
	for(uint32_t j = scores.size() - 1; j >= 0; --j){
		//Walk through column from top to bottom
		for(uint32_t i = 0; i < scores[j].size(); ++i){
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
			if(skT[col2pos[i]] == skT[col2pos[j]] && abcs[skT[col2pos[j]]] < 2) continue;

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