#ifndef SKETCH_TEST_HPP
#define SKETCH_TEST_HPP

#include "Sketch.h"

using namespace std;

class BuildSketchTest : public ::testing::Test {

	protected:

		BuildSketchTest(){}

		//The sketch to be calculated;
		Sketch s;
};

class SmHshsFrstTest : public ::testing::Test {

	protected:

		SmHshsFrstTest(): a(make_pair(1, 0)), b(make_pair(0, 1)) {}

		//One element
		pair<uint32_t, uint64_t> a;
		//Another element
		pair<uint32_t, uint64_t> b;
};

class RemDuplHshsTest : public ::testing::Test {

	protected:

		RemDuplHshsTest(): e(make_pair(0, 0)) {
			s.push_back(e);
		}

		//An element
		pair<uint32_t, uint64_t> e;
		//A sketch
		Sketch s;
};

#endif