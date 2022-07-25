#include "Sketch.cpp"
#include "IO.cpp"
#include "Thomology.h"

int main(int argc, char **argv){
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
	string seqA, seqB;
	//The input sequences' sketches
	Sketch skA, skB;

	//Parse arguments
	if(!prsArgs(argc, argv, seqA, seqB, kmerLen, frac, comWght, uniWght, tThres)){
		//Display help message
		dsHlp();
		return 1;
	}

	//Calculate sketches
	skA = buildSketch(seqA);
	skB = buildSketch(seqB);

	return 0;
}