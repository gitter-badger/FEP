#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>

#include "ElasticAnalysis2D.h"
#include "MeshBuilder.h"

class ElasticAnalysisTest : public testing::Test
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

TEST_F(ElasticAnalysisTest, FunctionReturns) {
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 2.0);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	apf::changeMeshShape(mesh, apf::getLagrange(2));

	ElasticAnalysis2D* tmp = new ElasticAnalysis2D(mesh);

	EXPECT_EQ(0, tmp->setup());
	EXPECT_EQ(0, tmp->solve());
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = mesh->begin(mesh->getDimension());
	while((e = mesh->iterate(it))) {
		EXPECT_EQ(0, tmp->makeStiffnessContributor(e));
		EXPECT_EQ(0, tmp->makeForceContributor(e));
		EXPECT_EQ(0, tmp->makeConstraint(e));
	}
	EXPECT_EQ(0, tmp->recover());

	delete tmp;
}
