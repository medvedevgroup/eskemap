#ifndef MEASURES_HPP
#define MEASURES_HPP

#include "Sketch.h"

enum Measure{none, intersec};

//This function calculates the "set intersection" similarity score
const int32_t calcIntersecScore(const Sketch& skA, const Sketch& skB);

#endif