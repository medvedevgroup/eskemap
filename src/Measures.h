#ifndef MEASURES_HPP
#define MEASURES_HPP

#include "Sketch.h"

enum Measure{none, intersec, algnWoutOffs};

//This function calculates the "set intersection" similarity score
const int32_t calcIntersecScore(const PairSketch& skA, const PairSketch& skB);

//This function calculates the "maximum aligned hashes" similarity score
const int32_t calcAlgnHshsScore(const PairSketch& skA, const PairSketch& skB, const bool& consOffs);

#endif