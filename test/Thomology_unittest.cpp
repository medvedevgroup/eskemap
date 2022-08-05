#include <gtest/gtest.h>

#include "ThomologyTest.h"
#include "Thomology.cpp"

//Tests for function const vector<Thomology> findThoms(const Sketch&, const Sketch&, const uint32_t&, const uint32_t&, const int32_t&)//
//	1. A hash is (not) unique in the pattern sketch DONE
//	2. Text and pattern sketches do (not) share a hash DONE
//	3. A hash inside the text sketch appears (not) for the first time DONE
//	4. While iterating over a list in pos, we dealed with a position k that was (not) k_min DONE
//	5. When calculating a score, the base case was (not) used DONE
//	6. While calculating a score recursively, we could (not) match a hash DONE
//	7. We are (not) dealing with a substring whose first and last hash are the same and that hash only occurs once inside the pattern DONE
//	8. The list of maximum scores contains a score found in a row that is (not) larger than the current row DONE
//	9. The score of a maximum t-homology does (not) give rise to an update of the threshold to compare against DONE
//	10. A t-homology is (not) maximal DONE

//Tests for function findThoms under the following conditions
//	1. A hash is (not) unique in the pattern sketch
//	2. Text and pattern sketches do (not) share a hash
//	3. A hash inside the text sketch appears (not) for the first time
//	4. While iterating over a list in pos, we dealed with a position k that was (not) k_min
//	5. When calculating a score, the base case was (not) used
//	6. While calculating a score recursively, we could (not) match a hash
//	7. We are (not) dealing with a substring whose first and last hash are the same and that hash only occurs once inside the pattern
//	8. The list of maximum scores contains a score found in a row that is (not) larger than the current row
//	9. The score of a maximum t-homology does give rise to an update of the threshold to compare against
//	10. A t-homology is (not) maximal
TEST_F(FindThomsTest, uniH){
	r = findThoms({1, 2, 1}, {1, 2, 3, 2, 2, 1}, 1, 1, -1);

	ASSERT_EQ(r.size(), 3);
	EXPECT_EQ(get<0>(r.front()), 0);
	EXPECT_EQ(get<1>(r.front()), 5);
	EXPECT_EQ(get<2>(r.front()), 0);
	i = r.begin();
	++i;
	EXPECT_EQ(get<0>(*i), 4);
	EXPECT_EQ(get<1>(*i), 5);
	EXPECT_EQ(get<2>(*i), 1);
	EXPECT_EQ(get<0>(r.back()), 0);
	EXPECT_EQ(get<1>(r.back()), 1);
	EXPECT_EQ(get<2>(r.back()), 1);
}

//Tests for function findThoms under the following conditions
//	1. A hash is unique in the pattern sketch
//	2. Text and pattern sketches do (not) share a hash
//	3. A hash inside the text sketch appears (not) for the first time
//	4. While iterating over a list in pos, we dealed with a position k that was (not) k_min
//	5. When calculating a score, the base case was (not) used
//	6. While calculating a score recursively, we could (not) match a hash
//	7. We are (not) dealing with a substring whose first and last hash are the same and that hash only occurs once inside the pattern
//	8. The list of maximum scores contains a score found in a row that is (not) larger than the current row
//	9. The score of a maximum t-homology does (not) give rise to an update of the threshold to compare against
//	10. A t-homology is (not) maximal
TEST_F(FindThomsTest, irrMax){
	r = findThoms({1, 2, 3, 4}, {1, 5, 3, 4, 1, 2, 5, 5, 2}, 1, 1, -1);

	ASSERT_EQ(r.size(), 4);
	EXPECT_EQ(get<0>(r.front()), 0);
	EXPECT_EQ(get<1>(r.front()), 8);
	EXPECT_EQ(get<2>(r.front()), -1);
	i = r.begin();
	++i;
	EXPECT_EQ(get<0>(*i), 2);
	EXPECT_EQ(get<1>(*i), 8);
	EXPECT_EQ(get<2>(*i), 1);
	++i;
	EXPECT_EQ(get<0>(*i), 0);
	EXPECT_EQ(get<1>(*i), 5);
	EXPECT_EQ(get<2>(*i), 2);
	EXPECT_EQ(get<0>(r.back()), 2);
	EXPECT_EQ(get<1>(r.back()), 5);
	EXPECT_EQ(get<2>(r.back()), 4);
}
