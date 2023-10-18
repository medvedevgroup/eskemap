#ifndef MEASURES_HPP
#define MEASURES_HPP

#include "Sketch.h"

enum Measure{none, intersec, algnWoutOffs};

//This function calculates the "set intersection" similarity score
const int32_t calcIntersecScore(const PairSketch& skA, const PairSketch& skB);

//This function calculates the "maximum aligned hashes" similarity score
const int32_t calcAlgnHshsScore(const PairSketch& skA, const PairSketch& skB, const bool& consOffs);

//This function calculates the linear score
inline const float calcLinScore(const uint32_t& xmin, const uint32_t& pLen, const uint32_t& i, const uint32_t& j, const float& uw){//TODO: This function still needs to be tested!
	return xmin - w * (max(pLen, j - i + 1) - xmin);
}

#endif
