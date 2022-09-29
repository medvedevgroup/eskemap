#ifndef INDEX_TEST_HPP
#define INDEX_TEST_HPP

#include "Index.h"

class IndexTest : public ::testing::Test {

	protected:

		IndexTest() {}

		//First hash position pair
		pair<uint64_t, uint32_t> fp;
		//Second hash position pair
		pair<uint64_t, uint32_t> sp;
};

#endif