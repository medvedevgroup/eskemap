// #include "mmpriv.h"
#include "Sketch.cpp"
#include "IO.cpp"
#include "Thomology.cpp"
#include "Index.cpp"

// #include "options.c"
// #include "kthread.c"
// #include "misc.c"
// #include "kalloc.c"
// #include "bseq.c"
// #include "index.c"
// #include "map.c"

// mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
// {
// 	int64_t is_idx;
// 	mm_idx_reader_t *r;
// 	is_idx = mm_idx_is_idx(fn);
// 	if (is_idx < 0) return 0; // failed to open the index
// 	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
// 	r->is_idx = is_idx;
// 	if (opt) r->opt = *opt;
// 	else mm_idxopt_init(&r->opt);
// 	if (r->is_idx) {
// 		r->fp.idx = fopen(fn, "rb");
// 		r->idx_size = is_idx;
// 	} else r->fp.seq = mm_bseq_open(fn);
// 	if (fn_out) r->fp_out = fopen(fn_out, "wb");
// 	return r;
// }

//The FracMinHash ratio
double hFrac = HASH_RATIO;
unordered_map<uint64_t, char> bLstmers;

int main(int argc, char **argv){
	//Flag to save that scores are to be normalized
	bool normalize = NORM_FLAG_DEFAULT;
	//Flag to state if we are interested in nested results
	bool noNesting = NESTING_FLAG_DEFAULT;
	//The k-mer length
	uint32_t kmerLen = K;
	//Scoring weights
	uint32_t comWght = DEFAULT_WEIGHT;
	float uniWght = DEFAULT_WEIGHT;
	//The t-homology threshold
	float tThres = T;
	//Intercept and decent to interpolate thresholds
	float dec = 0;
	float inter = 0;
	//Input file names
	string pFile, tFile, bLstFl;
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
	vector<tuple<string, uint32_t, Sketch>> pSks;
	//An iterator to iterate over pattern sketches
	vector<tuple<string, uint32_t, Sketch>>::const_iterator p;

	//Parse arguments
	if(!prsArgs(argc, argv, pFile, tFile, kmerLen, hFrac, bLstFl, comWght, uniWght, tThres, normalize, dec, inter, noNesting)){//TODO: Tests for this function need to be adapted!
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
	bLstmers[40064324] = 1;
	bLstmers[8867545] = 1;
	string genome;
	// cout << "main: tFile: " << tFile << endl;
	// unordered_map<uint64_t, char> seenHashes;
	readFASTA(tFile, genome);
	Sketch tsk = buildSketch(genome, kmerLen, hFrac, bLstmers);
	// cout << "main: Length of unfiltered text sketch: " << tsk.size() << endl;
	int nHits;
	for(Sketch::const_iterator gi = tsk.begin(); gi != tsk.end(); ++gi){
		// if(!seenHashes.contains(*gi)){
		// 	seenHashes[*gi] = 1;
			const uint64_t *idx_p = mm_idx_get(tidx, *gi, &nHits);
		// 	// cout << nHits << endl;
		// 	if(nHits <= 10)	cout << *gi << endl;
		// }
		if(nHits > 0){
			//Iterate over all occurrences
		    for(uint32_t i = 0; i < nHits; ++i){
		    	cout << "Position of k-mer " << *gi << ":" << (((uint32_t)(*idx_p))>>1) << endl;
		        //Move to next occurrence
		        idx_p++;
	        }
		}
	}
	return 0;

	//Load high abundance k-mers
	bLstmers = readBlstKmers("highAbundKmersLrgr10.txt");
	//Open stream to read in patterns
	fStr.open(pFile);

	//Testing
	// tsk = buildSketch(genome, kmerLen, hFrac, bLstmers);
	// cout << "main: Length of filtered text sketch: " << tsk.size() << endl;
	// return 0;

	//Load pattern sequences in batches
	while(lPttnSks(fStr, kmerLen, hFrac, bLstmers, pSks) || !pSks.empty()){//TODO: Test for this function need to be adaptated!
		//Iterate over pattern sketches
		for(p = pSks.begin(); p != pSks.end(); ++p){
			//Only output pattern sequence name if there is more than one sequence
			if(pSks.size() > 1) cout << get<0>(*p) << endl;

			//Testing
			// cout << "Pattern sketch:" << endl;
			// for(Sketch::const_iterator k = get<2>(*p).begin(); k != get<2>(*p).end(); ++k) cout << *k << " ";
			// cout << endl;

			//Calculate an adapted threshold if we have the necessary informations
			if(dec != 0 && inter != 0) tThres = dec * get<1>(*p) + inter;

			//Testing
			// cout << "main: pattern length: " << get<1>(*p) << endl;
			cout << "main: tThres: " << tThres << endl;
			// cout << "main: noNesting flag is " << (noNesting ? "" : "not ") << "set" << endl; 

			//Find t-homologies and output them
			outputHoms(findThoms(get<2>(*p), tidx, comWght, uniWght, tThres, noNesting), normalize, get<2>(*p).size());//TODO: Tests for this function need to be adaptated!
		}

		//Remove processed pattern sketches
		pSks.clear();
	}

	return 0;
}