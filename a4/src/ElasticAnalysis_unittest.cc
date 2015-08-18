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

TEST_F(ElasticAnalysisTest, AppRunTest) {
	mesh_builder->build2DRectQuadMesh(mesh, 2, 1, 0.0, 0.0, 2.0, 1.0);
	EXPECT_TRUE(mesh != NULL);
	//apf::changeMeshShape(mesh, apf::getSerendipity());
	uint32_t polynomial_order = 1;
	apf::changeMeshShape(mesh, apf::getLagrange(polynomial_order));
	/*physical parameters*/

	double E, Nu;
	E = 1e8;
	Nu = 0.35;
	uint32_t integration_order = 4;
	bool reorder_flag = true;
	struct ElasticAnalysisInput input = {
			mesh,
			integration_order,
			polynomial_order,
			E,
			Nu,
			reorder_flag};

	ElasticAnalysis2D tmp(input);

	EXPECT_EQ(0, tmp.setup());
	EXPECT_EQ(0, tmp.solve());
	EXPECT_EQ(0, tmp.recover());
}
