#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <gmi_null.h>


#include "MeshBuilder.h"
#include "MeshAdjReorder.h"

class ShapeFunctionOrder : public testing::Test
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;

	virtual void SetUp() {
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

TEST_F(ShapeFunctionOrder, Quadrilaterals) {
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 2.0);
	EXPECT_TRUE(mesh != NULL);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::changeMeshShape(mesh, apf::getLagrange(2));
	apf::Numbering* nodeNums = apf::createNumbering(mesh, "nodeNums", mesh->getShape(),1);
    apf::Numbering* faceNums = apf::createNumbering(mesh, "faceNums", apf::getConstant(mesh->getDimension()), 1);

	adjReorder(mesh, mesh->getShape(), 1, nodeNums, faceNums);

	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = mesh->begin(mesh->getDimension());
	while((e = mesh->iterate(it))) {
		uint32_t nnodes = apf::countElementNodes(mesh->getShape(), mesh->getType(e));
		apf::NewArray< int > node_mapping(nnodes);
		apf::getElementNumbers(nodeNums, e, node_mapping);
		for(uint32_t ii = 0; ii < nnodes; ++ii){
			std::cout << "Node " << ii << ": " << node_mapping[ii] << std::endl;
		}
		std::cout << "Mesh element processed" << std::endl;
	}
}
