#include <gtest/gtest.h>

#include "SketchTest.h"
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

//Tests for function uint64_t calcKmerNb(const string&)//
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
	EXPECT_EQ(getHash(27756, 262143), 197108);
}

//Tests for function PairSketch buildSketch(const string&)//
//	1. The input sequence length is smaller/equal/larger than/to k DONE
//	2. A hash value is (not) too large to be put into the sketch DONE

//Tests the function buildSketch under the following conditions
//	1. The input sequence length is smaller than k
TEST_F(BuildSketchTest, shrtSeq){
	s = buildSketch("a");
	EXPECT_TRUE(s.empty());
}

//Tests the function buildSketch under the following conditions
//	1. The input sequence length is equal to k
//	2. A hash value is too large to be put into the sketch
TEST_F(BuildSketchTest, lrgHsh){
	s = buildSketch("ACGTACGTA");
	EXPECT_TRUE(s.empty());
}

//Tests the function buildSketch under the following conditions
//	1. The input sequence length is larger than k
//	2. A hash value is (not) too large to be put into the sketch
TEST_F(BuildSketchTest, lngSeq){
	s = buildSketch("ACGTACGTACG");
	ASSERT_FALSE(s.empty());
	EXPECT_EQ(s.size(), 1);
	EXPECT_EQ(s.front().first, 1);
	EXPECT_EQ(s.front().second, getHash(calcKmerNb("CGTACGTAC"), pow(4, K) - 1));
}

//Tests for function inline const bool smHshsFrst(const pair<uint32_t, uint64_t>&, const pair<uint32_t, uint64_t>&)//
//	1. First element is (not) larger DONE

//Tests the function smHshsFrst under the following conditions
//	1. First element is larger
TEST_F(SmHshsFrstTest, FrstSml){
	EXPECT_TRUE(smHshsFrst(a, b));
}

//Tests the function smHshsFrst under the following conditions
//	1. First element is not larger
TEST_F(SmHshsFrstTest, FrstLrg){
	EXPECT_FALSE(smHshsFrst(b, a));
}

//Tests for function void remDuplHshs(PairSketch&)//
//	1. Sketch does (not) consist of more than 1 element DONE
//	2. An element can(not) be erased DONE

//Tests for function remDuplHshs under the following conditions
//	1. Sketch does consist of more than 1 element
//	2. An element can(not) be erased
TEST_F(RemDuplHshsTest, multElms){
	e = make_pair(1, 1);
	s.push_back(e);
	e = make_pair(2, 0);
	s.push_back(e);
	remDuplHshs(s);
	EXPECT_EQ(s.size(), 2);
	EXPECT_EQ(s.front().first, 0);
	EXPECT_EQ(s.front().second, 0);
	EXPECT_EQ(s.back().first, 1);
	EXPECT_EQ(s.back().second, 1);
}

//Tests for function remDuplHshs under the following conditions
//	1. Sketch does not consist of more than 1 element
TEST_F(RemDuplHshsTest, snglElm){
	remDuplHshs(s);
	EXPECT_EQ(s.size(), 1);
	EXPECT_EQ(s.front().first, 0);
	EXPECT_EQ(s.front().second, 0);
}

//Tests for function const Sketch buildSketch(const string&, const uint32_t&, const double&, const unordered_map<uint64_t, char>&)//
//	1. Sequence is (not) smaller than k DONE
//	2. A hash is (not) larger than the threshold DONE
//	3. A hash is (not) on the black list DONE

//Tests for function buildSketch under the following conditions
//	1. Sequence is smaller than k
TEST_F(BuildSketch1Test, shrtSeq){
	s = buildSketch("A", k, r, b);

	EXPECT_TRUE(s.empty());
}

//Tests for function buildSketch under the following conditions
//	1. Sequence is not smaller than k
//	2. A hash is (not) larger than the threshold
//	3. A hash is not on the black list
TEST_F(BuildSketch1Test, lngSeq){
	s = buildSketch("ACGTACGTACG", k, r, b);

	ASSERT_FALSE(s.empty());
	EXPECT_EQ(s.size(), 1);
	EXPECT_EQ(s.front(), getHash(calcKmerNb("CGTACGTAC"), (2 << (2 * K) - 1) - 1));
}

//Tests for function buildSketch under the following conditions
//	1. Sequence is not smaller than k
//	2. A hash is (not) larger than the threshold
//	3. A hash is (not) on the black list
TEST_F(BuildSketch1Test, blKm){
	b[13466] = 1;

	s = buildSketch("ACGTACGTACG", k, r, b);

	ASSERT_TRUE(s.empty());
}
