#include <gtest/gtest.h>
#include "../src/hpmath.h"
#include "../src/parallel.h"
#include <vector>

using namespace std;

// Tests the log(1+x) high precision function
TEST(LogOnePlusXTest, HighPrecisionTest) {
	EXPECT_DOUBLE_EQ(0, logopx(0));
}

// Tests the log(1+x) high precision function
TEST(, HighPrecisionTest) {
	EXPECT_DOUBLE_EQ(1, ompowopxn(-0.000495, 6728));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}