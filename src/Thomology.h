#ifndef THOMOLOGY_HPP
#define THOMOLOGY_HPP

#include <tuple>

#include "../minimap2/minimap.h"

#define T 0
#define DEFAULT_WEIGHT 1

using Thomology = tuple<uint32_t, uint32_t, int32_t>;

//This function finds all t-homologies of a text with respect to some pattern using dynamic programming
const vector<Thomology> findThoms(const Sketch& skP, const mm_idx_t *tidx, const uint32_t& cw, 
	const float& uw, const float& t, const bool& noNesting);

#endif
