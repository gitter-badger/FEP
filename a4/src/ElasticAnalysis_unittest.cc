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
#include "GeometryMappings.h"

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
	/*currently unused*/
	GeometryMappings* geo_map = new GeometryMappings();
	struct ElasticAnalysisInput input = {
			mesh,
			geo_map,
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

TEST_F(ElasticAnalysisTest, PlaneStrainComputation) {
	double E = 8e8;
	/*impossible for Nu to go over 0.5*/
	double Nu = 0.35;

	/*manually compute D using different method*/
	double mult_factor = E /((1+Nu)*(1- 2* Nu));

	double elm1 = mult_factor * (1 - Nu);
	double elm2 = mult_factor * (Nu);
	double elm3 = mult_factor * (1 - 2 *Nu) / 2;

	bool use_plane_stress = false;
	apf::Matrix< 3,3 > result_D = buildD(E, Nu, use_plane_stress);

	EXPECT_FLOAT_EQ(elm1, result_D[0][0]);
	EXPECT_FLOAT_EQ(elm2, result_D[0][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[0][2]);
	EXPECT_FLOAT_EQ(elm2, result_D[1][0]);
	EXPECT_FLOAT_EQ(elm1, result_D[1][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[1][2]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][0]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][1]);
	EXPECT_FLOAT_EQ(elm3, result_D[2][2]);
}

TEST_F(ElasticAnalysisTest, PlaneStressComputation) {
	double E = 8e8;
	/*impossible for Nu to go over 0.5*/
	double Nu = 0.35;

	/*manually compute D using different method*/
	double mult_factor = E /(1- Nu * Nu);

	double elm1 = mult_factor ;
	double elm2 = mult_factor * (Nu);
	double elm3 = mult_factor * (1 - Nu) / 2;

	bool use_plane_stress = true;
	apf::Matrix< 3,3 > result_D = buildD(E, Nu, use_plane_stress);

	EXPECT_FLOAT_EQ(elm1, result_D[0][0]);
	EXPECT_FLOAT_EQ(elm2, result_D[0][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[0][2]);
	EXPECT_FLOAT_EQ(elm2, result_D[1][0]);
	EXPECT_FLOAT_EQ(elm1, result_D[1][1]);
	EXPECT_FLOAT_EQ(0.0, result_D[1][2]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][0]);
	EXPECT_FLOAT_EQ(0.0, result_D[2][1]);
	EXPECT_FLOAT_EQ(elm3, result_D[2][2]);
}
