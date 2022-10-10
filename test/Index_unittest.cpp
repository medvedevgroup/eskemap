#include <gtest/gtest.h>

#include "IndexTest.h"

//Tests for function inline const bool smHshnPos(const pair<uint64_t, uint32_t>&, const pair<uint64_t, uint32_t>&)//
//	1. Hashes are equal and first pair's position is (not) smaller DONE
//	2. Hashes are not equal and first pair's position is (not) smaller DONE

//Tests the function smHshnPos under the following conditions
//	1. Hashes are equal and first pair's position is smaller
TEST_F(SmHshnPosTest, hefs){
	fp = make_pair(42, 20);
	sp = make_pair(42, 21);

	EXPECT_TRUE(smHshnPos(fp, sp));
}

//Tests the function smHshnPos under the following conditions
//	1. Hashes are equal and first pair's position is not smaller
TEST_F(SmHshnPosTest, hefl){
	fp = make_pair(42, 21);
	sp = make_pair(42, 20);

	EXPECT_FALSE(smHshnPos(fp, sp));
}

//Tests the function smHshnPos under the following conditions
//	2. Hashes are not equal and first pair's position is smaller
TEST_F(SmHshnPosTest, hufs){
	fp = make_pair(41, 20);
	sp = make_pair(42, 21);

	EXPECT_TRUE(smHshnPos(fp, sp));
}

//Tests the function smHshnPos under the following conditions
//	2. Hashes are not equal and first pair's position is not smaller
TEST_F(SmHshnPosTest, hufu1){
	fp = make_pair(41, 21);
	sp = make_pair(42, 21);

	EXPECT_TRUE(smHshnPos(fp, sp));
}

//Tests for function inline const bool smPos(const pair<uint64_t, uint32_t>&, const pair<uint64_t, uint32_t>&)//
//	1. First pair's position is (not) smaller 0/0

//Tests for function smPos under the following conditions
//	1. First pair's position is smaller
TEST(SmPosTest, fstSml){
	EXPECT_TRUE(smPos(make_pair(42, 23), make_pair(42, 24)));
}

//Tests for function smPos under the following conditions
//	1. First pair's position is not smaller
TEST(SmPosTest, eq){
	EXPECT_FALSE(smPos(make_pair(42, 23), make_pair(42, 23)));
}

//Tests for function const vector<pair<uint64_t, uint32_t>> genL(const Sketch&, const mm_idx_t*, const uint32_t&)//
//	1. The pattern sketch is (not) empty DONE
//	2. A hash can(not) be found inside the index DONE
//	3. A shared hash occurs (not) more than once inside the text sketch DONE
//	4. The order of shared hashes between pattern and text sketch is (not) different DONE

//Tests the function genL under the following conditions
//	1. The pattern sketch is empty
TEST_F(GenLtest, noP){
	L = genL(p, NULL, iopt.k);

	EXPECT_TRUE(L.empty());
}

//Tests the function genL under the following conditions
//	1. The pattern sketch is not empty
//	2. A hash can(not) be found inside the index
//	3. A shared hash occurs (not) more than once inside the text sketch
//	4. The order of shared hashes between pattern and text sketch is different
TEST_F(GenLtest, uH){
	L = genL({197108, 19367, 241811}, idx, iopt.k);

	ASSERT_FALSE(L.empty());
	ASSERT_EQ(L.size(), 3);
	EXPECT_EQ(L[0].first, 197108);
	EXPECT_EQ(L[0].second, 0);
	EXPECT_EQ(L[1].first, 241811);
	EXPECT_EQ(L[1].second, 3);
	EXPECT_EQ(L[2].first, 197108);
	EXPECT_EQ(L[2].second, 4);
}

//Tests the function genL under the following conditions
//	1. The pattern sketch is not empty
//	2. A hash can(not) be found inside the index
//	3. A shared hash occurs not more than once inside the text sketch
//	4. The order of shared hashes between pattern and text sketch is not different
TEST_F(GenLtest, eqOrd){
	L = genL({13466, 19367, 241811}, idx, iopt.k);

	ASSERT_FALSE(L.empty());
	ASSERT_EQ(L.size(), 2);
	EXPECT_EQ(L[0].first, 13466);
	EXPECT_EQ(L[0].second, 1);
	EXPECT_EQ(L[1].first, 241811);
	EXPECT_EQ(L[1].second, 3);
}
