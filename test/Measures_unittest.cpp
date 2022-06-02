#include <gtest/gtest.h>

#include "Measures.cpp"

//Tests for function const int32_t calcIntersecScore(const Sketch&, const Sketch&)//
//	1. The two sketches do (not) share an element 0/0

//Tests the function calcIntersecScore under the following conditions
//	1. The two sketches do (not) share an element
TEST(CalcIntersecScoreTest, SmplTst){
	Sketch a {make_pair(0, 43), make_pair(1, 44), make_pair(2, 42)};
	Sketch b {make_pair(0, 40), make_pair(1, 42), make_pair(2, 40)};

	EXPECT_EQ(calcIntersecScore(a, b), -2);
	EXPECT_EQ(a.size(), 3);
	Sketch::iterator i = a.begin();
	EXPECT_EQ(i->first, 2);
	EXPECT_EQ(i->second, 42);
	++i;
	EXPECT_EQ(i->first, 0);
	EXPECT_EQ(i->second, 43);
	++i;
	EXPECT_EQ(i->first, 1);
	EXPECT_EQ(i->second, 44);
	EXPECT_EQ(b.size(), 2);
	EXPECT_EQ(b.front().first, 0);
	EXPECT_EQ(b.front().second, 40);
	EXPECT_EQ(b.back().first, 1);
	EXPECT_EQ(b.back().second, 42);
}

//Tests for function const int32_t calcAlgnHshsScore(const Sketch&, const Sketch&, const bool&)//
//	1. We do (not) get a positive match score DONE
//	2. The maximum comes from a (mis)matching diagonal DONE
//	3. The maximum comes from above the current cell/the preceding cell DONE

//Tests the function calcAlgnHshsScore under the following conditions
//	1. We do (not) get a positive match score
//	2. The maximum comes from a (mis)matching diagonal
//	3. The maximum comes from above the current cell/the preceding cell
TEST(CalcAlgnHshsScore, SmplTst){
	Sketch a {make_pair(0, 0), make_pair(1, 2), make_pair(2, 3), make_pair(3, 3)};
	Sketch b {make_pair(0, 1), make_pair(1, 2), make_pair(2, 3), make_pair(3, 3)};

	EXPECT_EQ(calcAlgnHshsScore(a, b, false), 2);
}
