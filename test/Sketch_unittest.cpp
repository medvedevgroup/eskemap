#include <gtest/gtest.h>

//#include "SketchTest.h"
#include "../src/Sketch.cpp"

//Tests for function inline uint64_t getBaseNb(const char&)//
//	1. The character is A/C/G/T/not recognized DONE

//Tests the function getBaseNb under the following conditions
//	1. The character is A
TEST(GetBaseNbTest, A){
	EXPECT_EQ(getBaseNb('A'), 0);
}

//Tests the function getBaseNb under the following conditions
//	1. The character is C
TEST(GetBaseNbTest, C){
	EXPECT_EQ(getBaseNb('C'), 1);
}

//Tests the function getBaseNb under the following conditions
//	1. The character is G
TEST(GetBaseNbTest, G){
	EXPECT_EQ(getBaseNb('G'), 2);
}

//Tests the function getBaseNb under the following conditions
//	1. The character is T
TEST(GetBaseNbTest, T){
	EXPECT_EQ(getBaseNb('T'), 3);
}

//Tests the function getBaseNb under the following conditions
//	1. The character is a
TEST(GetBaseNbTest, a){
	EXPECT_EQ(getBaseNb('a'), 0);
}

//Tests for function uint64_t calcKmerNb(const std::string&)//
//	1. Simple test DONE


//Tests the function calcKmerNb under the following conditions
//	1. Simple test
TEST(CalcKmerNbTest, SmTst){
	EXPECT_EQ(calcKmerNb("ACGTACGTA"), 27756);
}

//Tests for function static inline uint64_t getHash(uint64_t, uint64_t)//
//	1. Simple test DONE

//Tests the function getHash under the following conditions
//	Simple test
TEST(GetHashTest, SmTst){
	EXPECT_EQ(getHash(27756, 262143), 42);
}