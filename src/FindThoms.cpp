#include "Sketch.cpp"
#include "IO.cpp"
#include "Thomology.h"

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
	double frac = HASH_RATIO;
	//Input file names
	string pFile, tFile;
	//An input sequence
	string seq;
	//The input sequences' sketches
	Sketch skP, skT;

	//Parse arguments
	if(!prsArgs(argc, argv, pFile, tFile, kmerLen, frac, comWght, uniWght, tThres, normalize)){
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
	skP = buildSketch(seq, kmerLen, frac);
	//We do not need the pattern sequence anymore
	seq.clear();

	//Try to load text sequence
	if(!readFASTA(tFile, seq)){
		cerr << "ERROR: Text sequence file could not be read" << endl;
		return -1;
	}

	//Calculate text sketch
	skT = buildSketch(seq, kmerLen, frac);
	//Find t-homologies and output them
	outputHoms(findThoms(skP, skT, comWght, uniWght, tThres, normalize)); //TODO: Implement this functions!

	return 0;
}