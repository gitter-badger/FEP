#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>

#include "ElasticAnalysis2D.h"

class ElasticAnalysisTest : public testing::Test
{
protected:
	apf::Mesh2* mesh;

	virtual void SetUp() {
		/*use mesh builder to abstract mesh creations*/
		mesh = NULL;
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}	
	}

};

TEST_F(ElasticAnalysisTest, FunctionReturns) {
	ElasticAnalysis2D* tmp = new ElasticAnalysis2D();

	EXPECT_EQ(0, tmp->setup());
	EXPECT_EQ(0, tmp->solve());
	EXPECT_EQ(0, tmp->makeStiffnessContributor(NULL));
	EXPECT_EQ(0, tmp->makeForceContributor(NULL));
	EXPECT_EQ(0, tmp->makeConstraint(NULL));
	EXPECT_EQ(0, tmp->recover());

	delete tmp;
}
