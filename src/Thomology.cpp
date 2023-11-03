#include "Thomology.h"
#include "Index.h"

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
void findThoms(const Sketch& skP, const mm_idx_t *tidx, const uint32_t& cw, const float& uw, const float& t, const bool& noNesting, 
	const bool& normalize){
	//A flag indicating whether the current candidate mapping is right reasonable
	bool isRr = true;
	//Some counter variables
	uint32_t i, occ, occPosJ;//, imin;
	//Index determining the position in xmin belonging to the last right reasonable candidate mapping
	uint32_t lstRrsnbl;
	//Another counter variable which may become negative
	int32_t j;
	//The maximum threshold to compare against
	float maxThres;
	//Current score
	float currScr;
	//A vector of x_min values of candidate mappings ending with the same k-mer
	vector<uint32_t> xmin;
	//An iterator to iterate over xmin's elements
	vector<uint32_t>::const_iterator vIt;
	vector<Thomology> res;
	
	// vector<Thomology>::const_reverse_iterator ri;

	//A sketch interator
	Sketch::const_iterator fSkIt;
	//Hash position array of all hashes and their positions in the text sketch which also appear inside the pattern
	vector<pair<uint64_t, uint32_t>> L;
	//An iterator to iterate over L
	vector<pair<uint64_t, uint32_t>>::const_iterator fLit;
	//Iterators to reversely iterate over L
	vector<pair<uint64_t, uint32_t>>::const_reverse_iterator lKmIt, fKmIt, jPosKmIt;
	//A list to store scores of final candidate mappings along with their start position
	list<pair<uint32_t, float>> maxScores;
	//An iterator for the maximum scores list
	list<pair<uint32_t, float>>::iterator li;
	//Initialize occp
	unordered_map<uint64_t, uint32_t> occp(skP.size());
	//Initialize hloc (again, we assume there are no duplicates in t)
	unordered_map<uint64_t, uint32_t> hloc(skP.size());

	//We do not support the exclusion of nested results currently
	if(noNesting){
		cerr << "ERROR: No nesting is currently not supported!" << endl;
		exit(-1);
	}

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
	//Reserve enough space for x_min vector
	xmin.reserve(L.size());

	//Base case: Deal with candidate mappings ending with last k-mer in L//

	//Set j to index of last k-mer in L
	j = L.size() - 1;
	//The smallest candidate mapping starting at a certain position is always right reasonable
	lstRrsnbl = j;
	//Set the threshold that a maximal candidate mapping has to exceed (-1, because once we found the first maximal candidate 
	//mapping, every other nested one needs to have a larger score)
	maxThres = t - 1;
	//Set iterator to iterate over last k-mer of each candidate mapping
	lKmIt = L.rbegin();
	//Add first count to hloc
	hloc[lKmIt->first] = 1;
	//x_min is 1 for the candidate mapping only containing the last k-mer of L by definition
	xmin.push_back(1);

	//Iterate over start positions in reverse order
	for(fKmIt = L.rbegin() + 1; fKmIt != L.rend(); ++fKmIt){//, ++itCounter
		//Fill hloc once iterating over L for the first time
		if(hloc.contains(fKmIt->first)){
			++hloc[fKmIt->first];
		} else{
			hloc[fKmIt->first] = 1;
		}

		//Check if candidate mapping is left reasonable
		if(occp[fKmIt->first] >= hloc[fKmIt->first]){
			xmin.insert(xmin.begin(), xmin.front() + 1);
		//Check if right reasonability no longer holds
		} else{
			xmin.insert(xmin.begin(), xmin.front());

			if(fKmIt->first == lKmIt->first) isRr = false;
		} 

		//Decrease right reasonability index if necessary
		lstRrsnbl -= (isRr ? 1 : 0);
	}

	//Find final mappings in last column
	for(i = 0; i < L.size(); ++i){
		//Calculate linear score for current mapping
		currScr = calcLinScore(xmin[i], skP.size(), L[i].second, L[j].second, uw);

		//Check if mapping is final
		if(currScr > maxThres && lstRrsnbl <= i){
			//A candidate mapping consisting of only one k-mer is always left reasonable; otherwise check if not left 
			//reasonable and continue
			if(i < L.size() - 1 && xmin[i] != xmin[i + 1] + 1) continue;

			//Add mapping to result list
			res.push_back(make_tuple(L[i].second, L[j].second, currScr));
			//Keep score for comparison to other scores of remaining candidate mappings starting with the same k-mer or later
			maxScores.push_back(make_pair(i, currScr));
			//Update score threshold 
			maxThres = currScr;
		}
	}

	//Iterate over all reasonable candidate mappings
	for(--j, jPosKmIt = L.rbegin() + 1; jPosKmIt != L.rend(); ++lKmIt, --j, ++jPosKmIt){
		//Initalize counter to count how many times the last k-mer of the previous candidate mapping has been seen while iterating 
		//over L from the beginning
		occ = 0;
		//Initalize counter to count how many times the last k-mer of the current candidate mapping has been seen while iterating 
		//over L from the beginning
		occPosJ = 0;
		//Reset flag
		isRr = false;
		//Update hloc
		--hloc[lKmIt->first];

		//Iterate over start positions in non-reverse order
		for(fLit = L.begin(), i = 0; i <= j; ++fLit, ++i){
			//Check right reasonability considering the previously last k-mer
			if(hloc[lKmIt->first] + 1 - occ <= occp[lKmIt->first]){
				//Decrement x_min
				xmin[i] -= 1;
			}

			//Check right reasonability considering the currently last k-mer
			if(hloc[jPosKmIt->first] - occPosJ <= occp[jPosKmIt->first]){
				//Save index in case this is the first right reasonable candidate mapping for the currently last k-mer we saw
				lstRrsnbl = (isRr ? lstRrsnbl : i);
				//Update flag
				isRr = true;
			}

			//occ changes for the next iteration if the first and the previously last k-mer of the candidate mapping have an 
			//identical sequence
			occ += (fLit->first == lKmIt->first ? 1 : 0);
			//occPosJ changes for the next iteration of the first and the currently last k-mer of the candidate mapping have an 
			//identical sequence
			occPosJ += (fLit->first == jPosKmIt->first ? 1 : 0);
		}
		
		//Iterate over start positions again
		for(fLit = L.begin(), i = 0, li = maxScores.begin(), maxThres = t - 1; i <= j; ++fLit, ++i){
			//Update threshold in case we have seen a relevant final mapping already
			if(li != maxScores.end() && li->first == i){
				maxThres = max(li->second, maxThres);
				//In memory of the "smooth iterator"
			}

			//Calculate linear score
			currScr = calcLinScore(xmin[i], skP.size(), L[i].second, L[j].second, uw);

			//Check if candidate mapping is final
			if(currScr > maxThres && i >= lstRrsnbl && (i == j || xmin[i] == xmin[i + 1] + 1)){
				//Add mapping to result list
				res.push_back(make_tuple(L[i].second, L[j].second, currScr));

				//Keep score for comparison to other scores of remaining candidate mappings starting with the same k-mer or later
				if(li != maxScores.end() && i == li->first){
					//Override existing maximum for candidate mapping starting with current k-mer
					*li = make_pair(i, currScr);
				} else{
					//Add a new entry
					maxScores.insert(li, make_pair(i, currScr));
				}

				//Update score threshold 
				maxThres = currScr;
			}

			//Move forward in maximum score list if needed
			if(li != maxScores.end() && li->first == i) ++li;
		}

		//Report results if we have already found a certain number
		if(res.size() > L.size()){
			outputHoms(res, normalize, skP.size());
			res.clear();
		}
	}

	outputHoms(res, normalize, skP.size());
}
