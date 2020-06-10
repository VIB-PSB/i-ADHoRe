#include <gtest/gtest.h>
#include "../src/DataSet.h"

TEST(indexToXYTest, 4x4Tests) {
	int x, y;

	DataSet::indexToXY(0, 4, x, y);
	EXPECT_DOUBLE_EQ(0, x);
	EXPECT_DOUBLE_EQ(0, y);

	DataSet::indexToXY(1, 4, x, y);
	EXPECT_DOUBLE_EQ(0, x);
	EXPECT_DOUBLE_EQ(1, y);

	DataSet::indexToXY(3, 4, x, y);
	EXPECT_DOUBLE_EQ(0, x);
	EXPECT_DOUBLE_EQ(3, y);

	DataSet::indexToXY(4, 4, x, y);
	EXPECT_DOUBLE_EQ(1, x);
	EXPECT_DOUBLE_EQ(1, y);

	DataSet::indexToXY(5, 4, x, y);
	EXPECT_DOUBLE_EQ(1, x);
	EXPECT_DOUBLE_EQ(2, y);

	DataSet::indexToXY(7, 4, x, y);
	EXPECT_DOUBLE_EQ(2, x);
	EXPECT_DOUBLE_EQ(2, y);

	DataSet::indexToXY(9, 4, x, y);
	EXPECT_DOUBLE_EQ(3, x);
	EXPECT_DOUBLE_EQ(3, y);
}

TEST(indexToXYTest, 2x2Tests) {
	int x, y;

	DataSet::indexToXY(0, 2, x, y);
	EXPECT_DOUBLE_EQ(0, x);
	EXPECT_DOUBLE_EQ(0, y);

	DataSet::indexToXY(1, 2, x, y);
	EXPECT_DOUBLE_EQ(0, x);
	EXPECT_DOUBLE_EQ(1, y);

	DataSet::indexToXY(2, 2, x, y);
	EXPECT_DOUBLE_EQ(1, x);
	EXPECT_DOUBLE_EQ(1, y);
}

TEST(indexToXYTest, 1x1Tests) {
	int x, y;

	DataSet::indexToXY(0, 1, x, y);
	EXPECT_DOUBLE_EQ(0, x);
	EXPECT_DOUBLE_EQ(0, y);
}
