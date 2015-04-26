#include <stdint.h>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfShape.h>

#include "MeshBuilder.h"

class MeshBuilderTest : public testing::Test
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

	virtual void SetUp() {
		/*use mesh builder to abstract mesh creations*/
		mesh_builder = new MeshBuilder();
		mesh = NULL;
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}	
		delete mesh_builder;
	}
};


TEST_F(MeshBuilderTest, PostitivePositive) {
	/*construct a mesh in first quadrant only*/
	mesh_builder->build2DRectQuadMesh(mesh, 10, 10, 1, 1, 20, 30);
	printf("Rest of test not implemented\n");
}

TEST_F(MeshBuilderTest, NegativeNegative) {
	/*create a mesh that starts in third and ends in third quadrant*/
	mesh_builder->build2DRectQuadMesh(mesh, 10, 10, -1, -1, -20, -30);
	printf("Rest of test not implemented\n");
}

TEST_F(MeshBuilderTest, PositiveNegative) {
	/*create a mesh that starts in first and ends in third quadrant*/
	mesh_builder->build2DRectQuadMesh(mesh, 10, 10, 1, 1, -20, -30);
	printf("Rest of test not implemented\n");
}

TEST_F(MeshBuilderTest, NegativePositive){
		/*create a mesh that starts in third and ends in first quadrant*/
	mesh_builder->build2DRectQuadMesh(mesh, 10, 10, -1, -1, 20, 30);
	printf("Rest of test not implemented\n");
}