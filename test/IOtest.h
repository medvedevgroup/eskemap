#ifndef IO_TEST_HPP
#define IO_TEST_HPP

#include "Measures.h"
#include "Thomology.h"
#include "IO.h"

using namespace std;

class PrsArgsTest : public ::testing::Test {

	protected:

		PrsArgsTest(): nbArgs(0), m(none) {}

		// void TearDown() override {
		// 	for(uint16_t i = 0; i < nbArgs; ++i) free(&argv[i]);
		// }

		//Measure to use
		Measure m;
		//Number of command line arguments
		int nbArgs;
		//Array with command line arguments
		char** argv;
		//First input sequence
		string sa;
		//Second input sequence
		string sb;
};

class PrsArgs1Test : public ::testing::Test {

	protected:

		PrsArgs1Test(): nbArgs(0), k(K), h(HASH_RATIO), c(DEFAULT_WEIGHT), u(DEFAULT_WEIGHT), t(T), n(NORM_FLAG_DEFAULT) {}

		//Normalization flag
		bool n;
		//k-mer length
		uint32_t k;
		//Common hash weight
		uint32_t c;
		//Unique hash weight
		uint32_t u;
		//Number of command line arguments
		int nbArgs;
		//t-homology threshold
		int32_t t;
		//Hash ratio
		double h;
		//Array with command line arguments
		char** argv;
		//Pattern sequence
		string p;
		//Text sequence
		string s;
};

#endif