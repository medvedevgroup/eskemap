#ifndef INDEX_TEST_HPP
#define INDEX_TEST_HPP

#include "Index.cpp"

double hFrac = 0.1;

class SmHshnPosTest : public ::testing::Test {

	protected:

		SmHshnPosTest() {}

		//First hash position pair
		pair<uint64_t, uint32_t> fp;
		//Second hash position pair
		pair<uint64_t, uint32_t> sp;
};

class GenLtest : public ::testing::Test {

	protected:

		GenLtest() {
			//Adjust FracMinHash ratio
			hFrac = 1.0;
			//Initialize index options
			mm_set_opt(0, &iopt, &mopt);
			//Adjust k if necessary
			iopt.k = K;
			//Open an index reader
			r = mm_idx_reader_open("testText.fasta", &iopt, INDEX_DEFAULT_DUMP_FILE);
			//Read index
			idx = mm_idx_reader_read(r, 1);
		}

		//An index option struct
		mm_idxopt_t iopt;
		//A mapping options struct
		mm_mapopt_t mopt;
		//An index reader
		mm_idx_reader_t *r;
		//A pointer to the index
		mm_idx_t* idx;
		//The pattern sketch
		Sketch p;
		//The L array
		vector<pair<uint64_t, uint32_t>> L;
};

#endif