#ifndef THOMOLOGY_TEST_HPP
#define THOMOLOGY_TEST_HPP

#include <vector>

#include "Sketch.h"
#include "Thomology.h"

class FindThomsTest : public ::testing::Test {

protected:

		FindThomsTest() {
			//Adjust FracMinHash ratio
			hFrac = 1.0;
			//Initialize index options
			mm_set_opt(0, &iopt, &mopt);
			//Adjust k if necessary
			iopt.k = K;
		}

		//An index option struct
		mm_idxopt_t iopt;
		//A mapping options struct
		mm_mapopt_t mopt;
		//An index reader
		mm_idx_reader_t *rd;
		//A pointer to the index
		mm_idx_t* idx;
		//The list of results
		vector<Thomology> r;
		//An iterator for the result list
		vector<Thomology>::const_iterator i;
};

#endif
