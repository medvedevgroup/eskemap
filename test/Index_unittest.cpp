#include <gtest/gtest.h>

#include "IndexTest.h"
#include "Index.cpp"

//Tests for function inline const bool smHshnPos(const pair<uint64_t, uint32_t>&, const pair<uint64_t, uint32_t>&)//
//	1. Hashes are equal and first pair's position is (not) smaller DONE
//	2. Hashes are not equal and first pair's position is (not) smaller DONE

//Tests the function smHshnPos under the following conditions
//	1. Hashes are equal and first pair's position is smaller
TEST_F(IndexTest, hefs){
	fp = make_pair(42, 20);
	sp = make_pair(42, 21);

	EXPECT_TRUE(smHshnPos(fp, sp));
}

//Tests the function smHshnPos under the following conditions
//	1. Hashes are equal and first pair's position is not smaller
TEST_F(IndexTest, hefl){
	fp = make_pair(42, 21);
	sp = make_pair(42, 20);

	EXPECT_FALSE(smHshnPos(fp, sp));
}

//Tests the function smHshnPos under the following conditions
//	2. Hashes are not equal and first pair's position is smaller
TEST_F(IndexTest, hufs){
	fp = make_pair(41, 20);
	sp = make_pair(42, 21);

	EXPECT_TRUE(smHshnPos(fp, sp));
}

//Tests the function smHshnPos under the following conditions
//	2. Hashes are not equal and first pair's position is not smaller
TEST_F(IndexTest, hufu1){
	fp = make_pair(41, 21);
	sp = make_pair(42, 21);

	EXPECT_TRUE(smHshnPos(fp, sp));
}
