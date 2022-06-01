#include <gtest/gtest.h>

#include "IOtest.h"
#include "IO.cpp"

//Tests for function const bool prsArgs(int&, char**, string&, string&, Measure&)//
//	1. Minimum number of parameters is (not) given DONE
//	2. First input sequence is (not) given DONE
//	3. Second input sequence is (not) given DONE
//	4. Similarity measure is (not) given 1/0
//	5. Help flag is (not) set DONE

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is not given
TEST_F(PrsArgsTest, fewPrms){
	nbArgs = 5;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[1] = strdup("-a");
	argv[2] = strdup("-b");
	argv[3] = strdup("ACGT");
	argv[4] = strdup("-i");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 5);
	EXPECT_EQ(sa, "");
	EXPECT_EQ(sb, "");
	EXPECT_EQ(none, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is set
TEST_F(PrsArgsTest, hlpFlg){
	nbArgs = 7;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[1] = strdup("-a");
	argv[2] = strdup("ACGT");
	argv[3] = strdup("-b");
	argv[4] = strdup("ACGT");
	argv[5] = strdup("-i");
	argv[6] = strdup("-h");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 7);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is not given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is not set
TEST_F(PrsArgsTest, noSeqA){
	nbArgs = 10;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[7] = strdup("-b");
	argv[8] = strdup("ACGT");
	argv[9] = strdup("-i");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 10);
	EXPECT_EQ(sa, "");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is not given
//	4. Similarity measure is given
//	5. Help flag is not set
TEST_F(PrsArgsTest, noSeqB){
	nbArgs = 13;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[10] = strdup("-a");
	argv[11] = strdup("ACGT");
	argv[12] = strdup("-i");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 13);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is not given
//	5. Help flag is not set
TEST_F(PrsArgsTest, noMes){
	nbArgs = 17;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[13] = strdup("-b");
	argv[14] = strdup("ACGT");
	argv[15] = strdup("-a");
	argv[16] = strdup("ACGT");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 17);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(none, m);
}
