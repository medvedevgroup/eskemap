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
	//The input sequences
	string seqP, seqT;
	//The input sequences' sketches
	Sketch skP, skT;

	//Parse arguments
	if(!prsArgs(argc, argv, seqP, seqT, kmerLen, frac, comWght, uniWght, tThres, normalize)){
		//Display help message
		dsHlp();
		return 1;
	}

	//Calculate sketches
	skP = buildSketch(seqP, kmerLen, frac);
	skT = buildSketch(seqT, kmerLen, frac);

	return 0;
}