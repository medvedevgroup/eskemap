#ifndef MEASURES_HPP
#define MEASURES_HPP

#include "Sketch.h"

enum Measure{none, intersec, algnWoutOffs};

//This function calculates the "set intersection" similarity score
const int32_t calcIntersecScore(const PairSketch& skA, const PairSketch& skB);

//This function calculates the "maximum aligned hashes" similarity score
const int32_t calcAlgnHshsScore(const PairSketch& skA, const PairSketch& skB, const bool& consOffs);

//This function calculates the linear score
inline const float calcLinScore(const int32_t& xmin, const int32_t& pLen, const int32_t& i, const int32_t& j, const float& w){//TODO: This function still needs to be tested!
	//Testing
	// if(i == 2501333 && j == 2503033) cout << "calcLinScore: xmin: " << xmin << " w: " << w << " pLen: " << pLen << " L[i]: " << i << " L[j]: " << j << endl;
	// cout << "calcLinScore: xmin: " << "L[j] - L[i]: " << j - i << endl;
	// cout << "calcLinScore: max(pLen, L[j] - L[i] + 1) - xmin: " << max(pLen, j - i + 1) - xmin << endl;

	return xmin - w * ((pLen + (j - i + 1)) - 2 * xmin);
}

#endif
