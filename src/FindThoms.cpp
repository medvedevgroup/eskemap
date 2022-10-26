#include "Sketch.cpp"
#include "IO.cpp"
#include "Thomology.cpp"
#include "Index.cpp"

//The FracMinHash ratio
double hFrac = HASH_RATIO;

int main(int argc, char **argv){
	//Flag to save that scores are to be normalized
	bool normalize = NORM_FLAG_DEFAULT;
	//The k-mer length
	uint32_t kmerLen = K;
	//Scoring weights
	uint32_t comWght = DEFAULT_WEIGHT;
	uint32_t uniWght = DEFAULT_WEIGHT;
	//The t-homology threshold
	int32_t tThres = T;
	//Input file names
	string pFile, tFile;
	//An input sequence
	string seq;
	//A file stream
	ifstream fStr;
	//A hash table to store black listed k-mers
	unordered_map<uint64_t, char> bLstmers;
	//An index option struct
	mm_idxopt_t iopt;
	//A mapping options struct
	mm_mapopt_t mopt;
	//An index reader
	mm_idx_reader_t *r;
	//A pointer to the index
	const mm_idx_t *tidx;
	//A vector of pattern sketches
	vector<pair<string, Sketch>> pSks;
	//An iterator to iterate over pattern sketches
	vector<pair<string, Sketch>>::const_iterator p;

	//Parse arguments
	if(!prsArgs(argc, argv, pFile, tFile, kmerLen, hFrac, comWght, uniWght, tThres, normalize)){
		//Display help message
		dsHlp();
		return 1;
	}

	//Set index options to default
	mm_set_opt(0, &iopt, &mopt);
	//Adjust k if necessary
	iopt.k = kmerLen;
	//Open an index reader //TODO: We do not allow yet to use a prebuilt index
	r = mm_idx_reader_open(tFile.c_str(), &iopt, INDEX_DEFAULT_DUMP_FILE);

	//Check if index could be opened successfully
	if(r == NULL){
		cerr << "ERROR: Text sequence file could not be read" << endl;
		return -1;
	}

	//Construct index
	if((tidx = mm_idx_reader_read(r, 1)) == 0){//TODO: Make use of multithreading here!
		cerr << "ERROR: Text index cannot be read" << endl;
		return -1;
	}

	//For simplicity we assume that an index always consists of only one part
	if(mm_idx_reader_read(r, 1) != 0){
		cerr << "ERROR: Text index consists of several parts! We cannot handle this yet" << endl;
		return -1; 
	}

	//Testing
	// string genome;
	// unordered_map<uint64_t, char> seenHashes;
	// readFASTA(tFile, genome);
	// Sketch tsk = buildSketch(genome, kmerLen, hFrac, bLstmers);
	// // cout << "main: Length of text sketch: " << tsk.size() << endl;
	// int nHits;
	// for(Sketch::const_iterator gi = tsk.begin(); gi != tsk.end(); ++gi){
	// 	if(!seenHashes.contains(*gi)){
	// 		seenHashes[*gi] = 1;
	// 		const uint64_t *idx_p = mm_idx_get(tidx, *gi, &nHits);
	// 		if(nHits > 100)	cout << *gi << endl;
	// 	}
	// }
	// return 0;

	//Load high abundance k-mers
	bLstmers = readBlstKmers("highAbundKmers.txt");//TODO: Implement this function!
	//Open stream to read in patterns
	fStr.open(pFile);

	//Load pattern sequences in batches
	while(lPttnSks(fStr, kmerLen, hFrac, bLstmers, pSks) || !pSks.empty()){//TODO: Implement this function!
		//Testing
		// cout << "main: Do we return?" << endl;
		// for(p = pSks.begin(); p != pSks.end(); ++p){
		// 	if(p->first == "S1_4" || p->first == "S1_5" || p->first == "S1_15"){
		// 		cout << p->first << endl;
		// 		int nHits;
		// 		for(Sketch::const_iterator hi = p->second.begin(); hi != p->second.end(); ++hi){
		// 			const uint64_t *idx_p = mm_idx_get(tidx, *hi, &nHits);
		// 			cout << *hi << " " << nHits << endl;
		// 		}
		// 	}
		// }
		// exit(0);

		//Iterate over pattern sketches
		for(p = pSks.begin(); p != pSks.end(); ++p){
			//Only output pattern sequence name if there is more than one sequence
			if(pSks.size() > 1) cout << p->first << endl;

			//Find t-homologies and output them
			outputHoms(findThoms(p->second, tidx, comWght, uniWght, tThres), normalize, p->second.size());
		}

		//Remove processed pattern sketches
		pSks.clear();

		//Testing
		// cout << "main: Cleared pattern sketch vector. Size: " << pSks.size() << " Empty: " << (pSks.empty() ? "Yes" : "No!") << endl;
	}

	return 0;
}