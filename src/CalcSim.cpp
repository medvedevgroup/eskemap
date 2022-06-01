#include "Sketch.cpp"
#include "Measures.cpp"
#include "IO.cpp"

int main(int argc, char **argv){
	//The measure to be used
	Measure msr = none;
	//The input sequences
	string seqA, seqB;
	//The input sequences' sketches
	Sketch skA, skB;

	//Parse arguments
	if(!prsArgs(argc, argv, seqA, seqB, msr)){
		//Display help message
		dspHlp();
		return 1;
	}

	//Calculate sketches
	skA = buildSketch(seqA);
	skB = buildSketch(seqB);

	//Calculate sketches' similarity
	cout << "Intersection similarity: " << calcIntersecScore(skA, skB) << endl;

	return 0;
}