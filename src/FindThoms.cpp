#include "Sketch.cpp"
#include "IO.cpp"
#include "Thomology.cpp"
#include "Index.cpp"

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
	//The FracMinHash ratio
	double hFrac = HASH_RATIO;
	//Input file names
	string pFile, tFile;
	//An input sequence
	string seq;
	//An index option struct
	mm_idxopt_t iopt;
	//A mapping options struct
	mm_mapopt_t mopt;
	//An index reader
	mm_idx_reader_t *r;
	//A pointer to the index
	const mm_idx_t *tidx;
	//The input sequences' sketches
	Sketch skP;
	//A vector storing hashes and their positions in the t which also occur in p
	vector<pair<uint64_t, uint32_t>> L;

	//Parse arguments
	if(!prsArgs(argc, argv, pFile, tFile, kmerLen, hFrac, comWght, uniWght, tThres, normalize)){
		//Display help message
		dsHlp();
		return 1;
	}

	//Try to load pattern sequence
	if(!readFASTA(pFile, seq)){
		cerr << "ERROR: Pattern sequence file could not be read" << endl;
		return -1;
	}

	//Calculate pattern sketch
	skP = buildSketch(seq, kmerLen, hFrac);
	//We do not need the pattern sequence anymore
	seq.clear();
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

	//Generate L
	L = genL(skP, tidx, kmerLen);//TODO: This function still needs to be tested!
	//Find t-homologies and output them
	outputHoms(findThoms(skP, L, comWght, uniWght, tThres), normalize, skP.size());//TODO: This function needs to be tested again!

	return 0;
}