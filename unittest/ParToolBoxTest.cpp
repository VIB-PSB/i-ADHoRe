#include <gtest/gtest.h>
#include <vector>
#include "../src/parallel.h"

using namespace std;

class ParallelTest : public ::testing::Test
{
protected:
	virtual void SetUp();

	// virtual void TearDown() {}

	vector<double> weights;
	vector<int> startPos1, startPos2, startPos3, startPos4;
};

void ParallelTest::SetUp()
{
	// simple partitioning setup
	for (int i = 1; i <= 10; i++)
		weights.push_back((double)i);

	ParToolBox::statPartitionWorkload(weights, startPos1, 4);

	// excess processes partitioning setup
	weights.clear();
	weights.push_back(1.0);
	weights.push_back(1.0);

	ParToolBox::statPartitionWorkload(weights, startPos2, 4);

	// empty weights setup
	weights.clear();

	ParToolBox::statPartitionWorkload(weights, startPos3, 4);

	// more difficult example
	weights.clear();
	weights.push_back(4);
	weights.push_back(4);
	weights.push_back(4);
	weights.push_back(7);
	weights.push_back(1);

	ParToolBox::statPartitionWorkload(weights, startPos4, 4);
}

TEST_F(ParallelTest, SimplePartitionTest) {
	EXPECT_EQ(4, startPos1.size());
	EXPECT_EQ(0, startPos1[0]);
	EXPECT_EQ(5, startPos1[1]);
	EXPECT_EQ(7, startPos1[2]);
	EXPECT_EQ(9, startPos1[3]);
}

TEST_F(ParallelTest, ExcessProcTest) {
	EXPECT_EQ(4, startPos2.size());
	EXPECT_EQ(0, startPos2[0]);
	EXPECT_EQ(1, startPos2[1]);
	EXPECT_EQ(2, startPos2[2]);
	EXPECT_EQ(2, startPos2[3]);
}

TEST_F(ParallelTest, EmptyWeightTest) {
	EXPECT_EQ(4, startPos3.size());
	EXPECT_EQ(0, startPos3[0]);
	EXPECT_EQ(0, startPos3[1]);
	EXPECT_EQ(0, startPos3[2]);
	EXPECT_EQ(0, startPos3[3]);
}

TEST_F(ParallelTest, PoisonedTest) {
	EXPECT_EQ(4, startPos4.size());
	EXPECT_EQ(0, startPos4[0]);
	EXPECT_EQ(1, startPos4[1]);
	EXPECT_EQ(2, startPos4[2]);
	EXPECT_EQ(3, startPos4[3]);
}
