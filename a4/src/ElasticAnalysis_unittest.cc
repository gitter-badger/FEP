#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>

#include "ElasticAnalysis2D.h"

class ElasticAnalysisTest : public testing::Test
{
protected:
	apf::Mesh2* mesh;

	virtual void SetUp() {
		gmi_register_null();
  		gmi_model* g = gmi_load(".null");
  		mesh = apf::makeEmptyMdsMesh(g, 2, false);
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}	
	}

};

TEST_F(ElasticAnalysisTest, FunctionReturns) {
	ElasticAnalysis2D* tmp = new ElasticAnalysis2D(mesh);

	EXPECT_EQ(0, tmp->setup());
	EXPECT_EQ(0, tmp->solve());
	EXPECT_EQ(0, tmp->makeStiffnessContributor(NULL));
	EXPECT_EQ(0, tmp->makeForceContributor(NULL));
	EXPECT_EQ(0, tmp->makeConstraint(NULL));
	EXPECT_EQ(0, tmp->recover());

	delete tmp;
}
