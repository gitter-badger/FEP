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

TEST_F(MatrixTest, Insertions) {
	BandedSymmetricMatrix* mat = new BandedSymmetricMatrix(10,10);

	std::cout << "test not implemented" << std::endl;
	std::cout << *mat << std::endl;

	delete mat;
}


