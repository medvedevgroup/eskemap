#ifndef IO_TEST_HPP
#define IO_TEST_HPP

#include "Measures.h"

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
		//Arrays with command line arguments
		char** argv;
		//First input sequence
		string sa;
		//Second input sequence
		string sb;
};

#endif