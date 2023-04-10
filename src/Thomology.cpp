#include "Thomology.h"
#include "Index.h"

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const mm_idx_t *tidx, const uint32_t& cw, 
	const float& uw, const float& t, const bool& noNesting){
	//Some counter variables
	uint32_t i, j, k, occ; //, maxI;
	//Length difference between the read's sketch and a homology interval
	int32_t lenDiff, minDiff;
	//The maximum threshold to compare against
	float maxThres;
	vector<Thomology> res;
	vector<Thomology>::const_reverse_iterator rResIt;
	//The score matrix
	vector<vector<float>> scores;
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
	list<tuple<uint32_t, float, int32_t>> maxScores;
	//An iterator for the maximum scores list
	list<tuple<uint32_t, float, int32_t>>::iterator li;
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

	//Generate L
	L = genL(occp, tidx);
	//Set position counter
	j = 0;

	//Fill score matrix
	for(fLit = L.begin(); fLit != L.end(); ++fLit){
		//Add new column in score matrix...
		scores.push_back(vector<float>());
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

	//Get a reverse iterator to iterate over L
	rLit = L.rbegin();
	//Get a counter for the column that we are at
	j = scores.size() - 1;

	// //Initialize last interesting row
	// maxI = scores.size();

	//Testing
	cout << "Find maximum t-homologies" << endl;

	//Find maximum t-homologies
	for(vector<vector<float>>::const_reverse_iterator colRit = scores.rbegin(); colRit != scores.rend(); ++colRit){
		//Reset row counter
		i = 0;
		//Start in the beginning of the list
		li = maxScores.begin();
		//Set the threshold that a maximum t-homology has to exceed
		maxThres = t - 1;
		//Get a forward iterator to iterate over L
		fLit = L.begin();
		minDiff = INT32_MAX;

		//Testing
		// cout << "j: " << j << endl;

		//Walk through column from top to bottom
		for(vector<float>::const_iterator rowIt = colRit->begin(); rowIt != colRit->end(); ++rowIt){
			//If the substring is a t-homology and its first and last hashes are identical this hash have lead to a positive contri-
			//bution to the score, otherwise this homology is not "reasonable"
			if(i != j && (scores[j - 1][i] < scores[j][i] - cw - uw + uw * (rLit->second - 
							L[j - 1].second - 1) || scores[j][i + 1] > scores[j][i] - cw - uw + uw * (L[i + 1].second - fLit->second
							 - 1))){
				++i;
				++fLit;
				continue;
			}

			//Testing
			if(j == 9230 && i == 6171) cout << "We are inside the relevant iteration" << endl;

			//Walk through the list to consider all relevant maximums seen so far
			while(li != maxScores.end()){
				//Maximums found in rows > i are irrelevant at this point
				if(i < get<0>(*li)) break;
			
				//Update maximum to compare with
				maxThres = max(maxThres, get<1>(*li));
				minDiff = get<2>(*li);
				//Walk on
				++li;
			}

			if(j == 9230 && i == 6171) cout << "Search for relevant maximums done" << endl;

			//We are not interested in any nested results
			if(noNesting && !maxScores.empty()){
				lenDiff = abs((int32_t) (L[j].second - L[i].second + 1 - skP.size()));

				//Testing
				// cout << "L[j].second: " << L[j].second << endl;
				// cout << "L[i].second: " << L[i].second << endl;
				// cout << "skP.size(): " << skP.size() << endl;
				// cout << "abs((int32_t) (L[j].second - L[i].second + 1 - skP.size())): " << abs((int32_t) (L[j].second - L[i].second + 1 - skP.size())) << endl;
				// cout << "(int32_t) (L[j].second - L[i].second + 1 - skP.size()) << 1 >> 1: " << ((int32_t) (L[j].second - L[i].second + 1 - skP.size()) << 1 >> 1) << endl;
				if(lenDiff < 0){
					cerr << "ERROR: lenDiff < 0. Something went wrong" << endl;
					exit(-1);
				}

				if(lenDiff > minDiff){
					++i;
					++fLit;
					continue;
				}

				// //There are no further interesting results in this column
				// if(i >= maxI) continue;

				// //Since we have no nested results we must not update our threshold
				// maxThres = t - 1;
			} //else{
			// 	//Walk through the list to consider all relevant maximums seen so far
			// 	while(li != maxScores.end()){
			// 		//Maximums found in rows > i are irrelevant at this point
			// 		if(i < get<0>(*li)) break;
			
			// 		//Update maximum to compare with
			// 		maxThres = max(maxThres, get<1>(*li));
			// 		//Walk on
			// 		++li;
			// 	}
			// }

			//Testing
			if(j == 9230 && i == 6171) cout << "Checked length of existing intervals" << endl;
			
			//Check if score is high enough
			if(*rowIt > maxThres){
				lenDiff = abs((int32_t) (L[j].second - L[i].second + 1 - skP.size()));
				//Add score to list with maximum scores
				maxScores.insert(li, make_tuple(i, *rowIt, lenDiff));

				//Testing
				if(j == 9230 && i == 6171) cout << "New maximum homology found; Deleting old stuff" << endl;

				//Delete homologies with a larger difference to the read if necessary
				if(noNesting){
					for(rResIt = res.rbegin(); rResIt != res.rend(); ++rResIt){
						if(get<0>(*rResIt) <= i){
							//Testing
							cout << "Deleting result?" << endl;

							res.erase(--(rResIt.base()));
						} else{
							break;
						}

						//Testing
						cout << "This is the last point where we get" << endl;
					}

					//Testing
					cout << "Leaving for" << endl;

					//Testing
					if(j == 9230 && i == 6171) cout << "Result list cleaned" << endl;

					li = maxScores.begin();

					while(li != maxScores.end() && !(get<0>(*li) == i && get<1>(*li) == *rowIt)){
						li = maxScores.erase(li);
					}
				}

				//Testing
				if(j == 9230 && i == 6171) cout << "Everything deleted" << endl;

				//Add t-homology to results
				res.push_back(make_tuple(i, j, *rowIt));
				//Update maximum to compare with
				maxThres = *rowIt;

				// //Update last interesting row (in case we do not want nested results)
				// maxI = i;

				minDiff = lenDiff;
			}

			++i;
			++fLit;
		}

		++rLit;
		--j;
	}

	//Testing
	cout << "Finished maximum search" << endl;

	//Map result coordinates to real positions inside the sketch
	for(vector<Thomology>::iterator ri = res.begin(); ri != res.end(); ++ri){
		//Testing
		cout << "Correcting coordinates" << endl;

		get<0>(*ri) = L[get<0>(*ri)].second;
		get<1>(*ri) = L[get<1>(*ri)].second;
	}

	//Testing
	cout << "But we do not get here" << endl;

	return res;
}
