#include <stdint.h>
#include <stdio.h>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfField.h>

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
	//apf::writeVtkFiles("outQuad", mesh);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::changeMeshShape(mesh, apf::getLagrange(1));
	//apf::writeVtkFiles("secondQuad", mesh);
	apf::MeshEntity* e;
	apf::MeshIterator* it = mesh->begin(2);
	//apf::Field * master_f = createFieldOn(mesh, "foo", apf::VECTOR); 
	apf::Field* master_f = createFieldOn(mesh, "foo", apf::SCALAR);
	apf::FieldShape* fs = mesh->getShape();

	while((e = mesh->iterate(it))) {
		//apf::MeshElement* mesh_elm = apf::getMeshEle
		//apf::Element* element_ptr = apf::createElement(master_f, e );
		apf::Vector3 param;
		double xyz[3] = {0, 0 , 0};
		param.fromArray(xyz);
		int order = 2;
		int rule = order * order;
		// for(int ii = 0; ii < rule; ++ ii) {
		// 	apf::getGaussPoint(, order, ii, param);
		// 	//std::cout << param << std::endl;
		// }
		std::cout << apf::countElementNodes(fs, mesh->getType(e)) << std::endl;
		//std::cout << master_f->countComponents() << std::endl;
		apf::EntityShape* es = fs->getEntityShape(mesh->getType(e));
		std::cout << "Order: " <<  fs->getOrder() << std::endl;
		std::cout << fs->countNodesOn(apf::Mesh::EDGE) << std::endl;
	}
}

TEST_F(RectMeshTest, Triangle) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	apf::writeVtkFiles("outTri", mesh);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("secondTri", mesh);
}
