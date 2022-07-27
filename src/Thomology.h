#ifndef THOMOLOGY_HPP
#define THOMOLOGY_HPP

#include <tuple>

#define T 0
#define DEFAULT_WEIGHT 1

using Thomology = tuple<uint32_t, uint32_t, int32_t>;

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const Sketch& skT, const uint32_t& cw, const uint32_t& uw, const int32_t& t, 
	const bool& norm);

#endif