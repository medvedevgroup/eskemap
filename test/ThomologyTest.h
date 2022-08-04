#ifndef THOMOLOGY_TEST_HPP
#define THOMOLOGY_TEST_HPP

#include <vector>

#include "Sketch.h"
#include "Thomology.h"

class FindThomsTest : public ::testing::Test {

protected:

		FindThomsTest() {}

		//The list of results
		vector<Thomology> r;
		//An iterator for the result list
		vector<Thomology>::const_iterator i;
};

#endif
