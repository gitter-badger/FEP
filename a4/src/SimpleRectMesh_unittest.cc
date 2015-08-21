#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfField.h>
#include <apfNumbering.h>
#include <apfMatrix.h>

#include "MeshBuilder.h"
#include "MeshAdjReorder.h"
#include "GeometryMappings.h"


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
	EXPECT_NO_THROW(mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 2.0));
	EXPECT_TRUE(mesh != NULL);
}

TEST_F(RectMeshTest, CheckBoundaryEdgeTags) {
	int X_ELMS = 3;
	int Y_ELMS = 3;
	EXPECT_NO_THROW(mesh_builder->build2DRectQuadMesh(mesh, X_ELMS, Y_ELMS, 0.0, 0.0, 2.0, 2.0));
	EXPECT_TRUE(mesh != NULL);
	apf::changeMeshShape(mesh, apf::getLagrange(1));
	apf::writeVtkFiles("testQuad", mesh);
	/*we do a bad thing and use undocumented numbering
	* that the mesh_builder creates to verify that numberings
	* are in the right place*/
	apf::MeshTag* edge_tag = mesh->findTag(EDGE_BC_TAG_NAME);
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = mesh->begin(1);
	while((e = mesh->iterate(it))) {
		apf::Downward down;
		int data;
		if(mesh->hasTag(e, edge_tag)) {
			/*find the vertices bounding that edge*/
			int n_verts = mesh->getDownward(e, 0, down);

		}
	}
	mesh->end(it);
}

TEST_F(RectMeshTest, CheckBoundaryVertexTags) {
	int X_ELMS = 3;
	int Y_ELMS = 3;
	EXPECT_NO_THROW(mesh_builder->build2DRectQuadMesh(mesh, X_ELMS, Y_ELMS, 0.0, 0.0, 2.0, 2.0));
	EXPECT_TRUE(mesh != NULL);
	apf::changeMeshShape(mesh, apf::getLagrange(1));
	apf::writeVtkFiles("testQuad", mesh);
	/*we do a bad thing and use undocumented numbering
	* that the mesh_builder creates to verify that numberings
	* are in the right place*/
	apf::MeshTag* vert_tag = mesh->findTag(VERT_BC_TAG_NAME);
	apf::Numbering* secret_numbering = mesh->findNumbering(SECRET_BUILDER_NUMBERING);
	EXPECT_TRUE(secret_numbering != NULL);
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = mesh->begin(0);
	while((e = mesh->iterate(it))) {
		int data;
		if(mesh->hasTag(e, vert_tag)) {
			mesh->getIntTag(e, vert_tag, &data);
			int node_label = apf::getNumber(secret_numbering, e, 0,0);
			/*we make sure that each corner of the rectangle is identifiable*/
			if( ((X_ELMS+1)*(Y_ELMS+1))== node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X1Y1) == data);
			} else if(1 == node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X0Y0) == data);
			} else if((X_ELMS * (Y_ELMS + 1) +1) == node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X0Y1) == data);
			} else if((X_ELMS * 2 + 1) == node_label) {
				EXPECT_TRUE(static_cast<int>(VERT_X1Y0) == data);
			}
		}
	}
	mesh->end(it);
}


TEST_F(RectMeshTest, Triangle) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	//apf::writeVtkFiles("outTri", mesh);
	apf::changeMeshShape(mesh, apf::getLagrange(2));
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("secondTri", mesh);
}
