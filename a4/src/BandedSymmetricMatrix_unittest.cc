#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include "BandedSymmetricMatrix.h"

class MatrixTest : public testing::Test
{
protected:

	virtual void SetUp(){

	}

	virtual void TearDown() {

	}
};

TEST_F(MatrixTest, AccessOperators) {
	BandedSymmetricMatrix mat(10,10);

	/*test for symmetry*/
	mat(0,9) = 27.0;
	EXPECT_EQ(27.0, mat(0,9));
	EXPECT_EQ(27.0, mat(9,0));

	/*test that non instantiated things off diagonal are zero,
	* we just quickly test some random indices in range*/
	EXPECT_EQ(0.0, mat(7,3));
	EXPECT_EQ(0.0, mat(2,3));
	EXPECT_EQ(0.0, mat(5,9));

}


