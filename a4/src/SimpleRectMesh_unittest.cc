#include <stdint.h>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfShape.h>

#include "MeshBuilder.h"


class RectMeshTest : public testing::Test
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

TEST_F(RectMeshTest, Rectangle) {
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0, 0, 2, 1);
	apf::writeVtkFiles("outQuad", mesh);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("secondQuad", mesh);
}

TEST_F(RectMeshTest, Triangle) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	apf::writeVtkFiles("outTri", mesh);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("secondTri", mesh);
}
