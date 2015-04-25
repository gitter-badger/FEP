#include <stdint.h>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>

#include "MeshBuilder.h"


class RectMeshTest : public testing::Test
{
protected:
	static const uint32_t x_units = 10;
	static const uint32_t y_uints = 10;

	virtual void SetUp() {
		printf("virtual void setup is called\n");
		apf::Mesh2* mesh;
		MeshBuilder mb = MeshBuilder();
		mb.build2DRectQuadMesh(mesh, 4, 5, 0, 0, 10, 20);
	}

	static int timesSeven(int n) {
		return n * 7;
	}
};

TEST_F(RectMeshTest, FirstTest) {
	EXPECT_EQ(7, timesSeven(3));
}

TEST_F(RectMeshTest, SecondTest) {
	EXPECT_EQ(14, timesSeven(2));
}
