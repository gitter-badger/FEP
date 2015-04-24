

#include <gtest/gtest.h>
#include <sys/types.h>

class RectMeshTest : public testing::Test
{
protected:
	uint32_t x_units = 10;
	uint32_t y_uints = 10;

	virtual void SetUp() {
		printf("virtual void setup is called")

	}

	static int timesSeven(int n) {
		return n * 7;
	}
};

TEST_F(RectMeshTest, FirstTest) {
	EXPECT_EQ(7, timesSeven(1));
}

TEST_F(RectMeshTest, SecondTest) {
	EXPECT_EQ(14, timesSeven(2))
}