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

	//Testing
	// cout << "findThoms: L.size(): " << L.size() << endl;

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

	//Testing
	// fKmIt = L.rbegin();
	// ++fKmIt;
	// cout << "findThoms: fKmIt->second (using rbegin() + ++): " << fKmIt->second << endl;
	// uint32_t itCounter = 0;

	//Iterate over start positions in reverse order
	for(fKmIt = L.rbegin() + 1; fKmIt != L.rend(); ++fKmIt){//, ++itCounter
		//Testing
		// cout << "findThoms: fKmIt->second (using rbegin() + 1): " << fKmIt->second << endl;
		// exit(0);

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

	//Testing
	// cout << "Filled first column of xmin" << endl;

	//Find final mappings in last column
	for(i = 0; i < L.size(); ++i){
		//Testing
		// cout << "i: " << i << endl;

		//Calculate linear score for current mapping
		currScr = calcLinScore(xmin[i], skP.size(), L[i].second, L[j].second, uw);

		//Testing
		// cout << "lstRrsnbl: " << lstRrsnbl << endl;
		// cout << "currScr: " << currScr << endl;
		// cout << "maxThres: " << maxThres << endl;

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

	//Testing
	// cout << "First column done" << endl;
	// cout << "Number of iterations to iterate over it: " << itCounter << endl;
	// cout << "res.size(): " << res.size() << endl;
	// cout << "maxThres: " << maxThres << endl;
	// cout << "xmin.size(): " << xmin.size() << endl;
	// for(uint32_t k = 0; k < 10; ++k){
	// 	cout << "xmin[" << k << "]: " << xmin[k] << endl;
	// }
	// // for(uint32_t k = 0; k < 10; ++k){
	// // 	cout << "L[" << k << "]: " << L[k].first << " " << L[k].second << endl;
	// // }
	// cout << "xmin[L.size()-2]: " << xmin[L.size()-2] << endl;
	// cout << "xmin[L.size()-1]: " << xmin[L.size()-1] << endl;
	// cout << "L[L.size()-2]: " << L[L.size()-2].first << " " << L[L.size()-2].second << endl;
	// cout << "L[L.size()-1]: " << L[L.size()-1].first << " " << L[L.size()-1].second << endl;
	// exit(0);

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

		//Testing
		// int32_t xminDecCount = xmin[i];

		//Iterate over start positions in non-reverse order
		for(fLit = L.begin(), i = 0; i <= j; ++fLit, ++i){
			//Testing
			// if((i == 5749 || i == 5748) && j == 9222){
			// 	cout << "i: " << i << endl;
			// 	cout << "hloc[jPosKmIt->first]: " << hloc[jPosKmIt->first] << endl;
			// 	cout << "occPosJ: " << occPosJ << endl;
			// 	cout << "occp[jPosKmIt->first]: " << occp[jPosKmIt->first] << endl;
			// }

			//Check right reasonability considering the previously last k-mer
			if(hloc[lKmIt->first] - occ <= occp[lKmIt->first]){
				//Decrement x_min
				xmin[i] -= 1;

				//Testing
				// --xminDecCount;
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

			//Testing
			// if(occ > hloc[lKmIt->first]){
			// 	cout << "occ > hloc[lKmIt->first]; this should never happen!" << endl;
			// 	cout << "occ: " << occ << endl;
			// 	cout << "hloc[" << lKmIt->first << "]: " << hloc[lKmIt->first] << endl;
			// 	cout << "i: " << i << endl;
			// 	cout << "j: " << j << endl;
			// 	// exit(0);
			// }
		}

		//Testing
		// if(xminDecCount < 0){
		// 	cout << "We decreased xmin too much" << endl;
		// 	cout << "i: " << i << " j: " << j << endl;
		// 	exit(0);
		// }

		//Update hloc
		--hloc[lKmIt->first];
		
		//Iterate over start positions again
		for(fLit = L.begin(), i = 0, li = maxScores.begin(), maxThres = t - 1; i <= j; ++fLit, ++i){
			//Testing
			// if(L[i].second == 3573864 && L[j].second == 3575170){
			// 	cout << "maxThres: " << maxThres << endl;
			// 	cout << "li->first: " << li->first << endl;
			// }

			//Update threshold in case we have seen a relevant final mapping already
			if(li != maxScores.end() && li->first == i){
				maxThres = max(li->second, maxThres);
				//In memory of the "smooth iterator"
			}

			//Calculate linear score
			currScr = calcLinScore(xmin[i], skP.size(), L[i].second, L[j].second, uw);

			//Testing
			// if(L[i].first == 1463247 && L[j].first == 1466043){
			// 	cout << "i: " << i << " j: " << j << endl;
			// 	cout << "currScr: " << currScr << endl;
			// 	cout << "maxThres: " << maxThres << endl;
			// 	cout << "lstRrsnbl: " << lstRrsnbl << endl;
			// 	cout << "xmin[i]: " << xmin[i] << endl;
			// 	cout << "xmin[i+1]: " << xmin[i+1] << endl;
			// 	// uint32_t lKmCnt = 0;
			// 	// for(uint32_t k = 6171; k < 9223; ++k){
			// 	// 	if(L[k].first == L[9222].first) ++lKmCnt;
			// 	// }
			// 	// cout << "Abundance of last k-mer in interval: " << lKmCnt << endl;
			// 	// cout << "Abundance of last k-mer in pattern: " << occp[L[9222].first] << endl;
			// 	// cout << "L[8189].first: " << L[8189].first << endl;
			// 	// cout << "L[9222].first: " << L[9222].first << endl;
			// }

			//Check if candidate mapping is final
			if(currScr > maxThres && i >= lstRrsnbl && (i == j || xmin[i] == xmin[i + 1] + 1)){
				//Testing
				// cout << "Final mapping found" << endl;
				// cout << "currScr: " << currScr << endl;
				// cout << "maxThres: " << maxThres << endl;
				// cout << "findThoms: i: " << i << " j: " << j << endl;
				// cout << "findThoms: xmin[i]: " << xmin[i] << endl;
				// exit(0);
				// if(L[i].second == 3573864 && L[j].second == 3575170){
				// 	cout << "currScr: " << currScr << endl;
				// 	cout << "maxThres: " << maxThres << endl;
				// 	cout << "L[i].first: " << L[i].first << endl;
				// 	cout << "i: " << i << " j: " << j << endl;
				// 	cout << "First hash in pattern sketch: " << skP.front() << endl;
				// 	cout << "skP.size(): " << skP.size() << endl;
				// 	exit(0);
				// }

				//Add mapping to result list
				res.push_back(make_tuple(L[i].second, L[j].second, currScr));

				//Keep score for comparison to other scores of remaining candidate mappings starting with the same k-mer or later
				if(i == li->first){
					//Override existing maximum for candidate mapping starting with current k-mer
					*li = make_pair(i, currScr);
				} else{
					//Testing
					// if(L[i].second == 3573864 && L[j].second == 3575171){
					// 	cout << "If we were here, we did not update but inserted a new maxium" << endl;
					// }

					//Add a new entry
					maxScores.insert(li, make_pair(i, currScr));
				}

				//Testing
				// if(L[i].second == 3573864 && L[j].second == 3575171){
				// 	cout << "Here we should have updated the score:" << endl;
				// 	cout << "currScr: " << currScr << endl;
				// 	cout << "li->second: " << li->second << endl;
				// }

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

		//Testing
		// cout << "Processing of column " << j << " done" << endl;
	}

	//Testing
	// cout << "Outputting results...res: " << endl;
	// cout << "res.size(): " << res.size() << endl;

	outputHoms(res, normalize, skP.size());

	//Testing
	// cout << "Read processed" << endl;
}
