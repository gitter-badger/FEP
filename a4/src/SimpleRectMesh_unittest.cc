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

#include "IntegrationHelper.h"
#include "TestClass.h"
#include "StiffnessContributor2D.h"


double foo_func(apf::Vector3 const& p) {
	return 0.0;
}

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
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 2.0);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::changeMeshShape(mesh, apf::getLagrange(1));
	apf::writeVtkFiles("before_secondQuad", mesh);

	apf::MeshEntity* e;
	apf::MeshIterator* it;	
	apf::FieldShape* fs = mesh->getShape();

	apf::Numbering* nodes_numbers = apf::numberOwnedNodes(mesh, "fatman", fs);

	//apf::Numbering* all_node_nums = apf::createNumbering(mesh, "allNodes", mesh->getShape(), 1);
	
    apf::Numbering* nodeNums = apf::createNumbering(mesh, "nodeNums", mesh->getShape(),1);
    apf::Numbering* faceNums = apf::createNumbering(mesh, "faceNums", apf::getConstant(mesh->getDimension()), 1);

	adjReorder(mesh, fs, 1, nodeNums, faceNums);
	// std::cout << nodeNums << std::endl;
	

	apf::Numbering* arb_nums = apf::createNumbering(mesh, "arb", mesh->getShape(), 2);

	apf::writeVtkFiles("secondQuad", mesh);

	int dim = mesh->getDimension();
	int batman = 0;	

	

	apf::Field* master_f = createField(mesh, "master_f", apf::VECTOR, mesh->getShape());
	apf::zeroField(master_f);

	apf::Matrix< 3,3 > D;


	StiffnessContributor2D integrate(master_f, D, 4);

	it = mesh->begin(2);
	while((e = mesh->iterate(it))) {
		//apf::MeshElement* mesh_elm = apf::getMeshElement()
		//apf::Element* element_ptr = apf::createElement(master_f, e );
		//apf::EntityShape::getValues
		apf::Vector3 param;
		double zeros[3] = {0, 0, 0};

		param.fromArray(zeros);
		int order = 4; //this is the order accuracy, not the polynomial degree

		apf::MeshElement* mesh_elm = apf::createMeshElement(mesh, e);
		integrate.process(mesh_elm);
		/*only print out the elements after processing*/
		std::cout << integrate.ke << std::endl;

		std::cout << "==============================" << std::endl;
		// apf::Element* field_elm = apf::createElement(master_f, mesh_elm);
		// std::cout << "Ndofs " << apf::countNodes(field_elm) << std::endl;
		// int num_components = master_f->countComponents();
		// std::cout << "num_components: " << num_components << std::endl;



		int num_int_points = apf::countIntPoints(mesh_elm, order);
		//apf::EntityShape* es = master_f->getEntityShape(mesh->getType(e));
		uint32_t num_nodes = apf::countElementNodes(fs, mesh->getType(e));

		apf::NewArray< int > node_mapping(num_nodes);

		apf::getElementNumbers(nodes_numbers, e, node_mapping);

		#if 0
		for(int ii = 0; ii < num_nodes; ++ii) {
			//std::cout << ii << " => " << node_mapping[ii] << std::endl;
		}

		for(int ii = 0; ii < num_int_points; ++ ii) {
			double running_sum = 0;

			apf::getIntPoint(mesh_elm, order, ii, param);
			//std::cout << std::setprecision(16);
			//std::cout << ii << ": " << param.x() << " " << param.y() << std::endl;

			apf::NewArray< apf::Vector3 > grads(num_nodes);
			es->getLocalGradients(mesh, e, param, grads);

			for(uint32_t jj = 0; jj < num_nodes; ++jj ) {

			}
		}
		std::cout << apf::countElementNodes(fs, mesh->getType(e)) << std::endl;
		//std::cout << master_f->countComponents() << std::endl;
		#endif
	}

}

TEST_F(RectMeshTest, Triangle) {
	mesh_builder->build2DRectTriMesh(mesh, 2, 1, 0, 0, 2, 1);
	apf::writeVtkFiles("outTri", mesh);
	apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::writeVtkFiles("secondTri", mesh);
}
