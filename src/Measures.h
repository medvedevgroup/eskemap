#ifndef MEASURES_HPP
#define MEASURES_HPP

#include "Sketch.h"

enum Measure{none, intersec, algnWoutOffs};

//This function calculates the "set intersection" similarity score
const int32_t calcIntersecScore(const Sketch& skA, const Sketch& skB);

//This function calculates the "maximum aligned hashes" similarity score
const int32_t calcAlgnHshsScore(const Sketch& skA, const Sketch& skB, const bool& consOffs);

#endif