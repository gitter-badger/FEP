#ifndef ELASTIC_ANALYSIS
#define ELASTIC_ANALYSIS 1

#include <apf.h>
#include <apfMesh.h>
#include <apfMatrix.h>
#include <apfNumbering.h>

#include "FEAnalysis.h"
#include "AlgebraicSystem.h"

struct ElasticAnalysisInput {
	apf::Mesh* m;
	uint32_t integration_order;
	uint32_t poly_order;
	double E;
	double Nu;
	bool reorder;
};

class ElasticAnalysis2D : FEAnalysis
{
public:
	ElasticAnalysis2D(struct ElasticAnalysisInput & in);
	virtual ~ElasticAnalysis2D();

	virtual uint32_t setup();
	virtual uint32_t solve();
	virtual uint32_t makeStiffnessContributor(apf::MeshEntity* e);
	virtual uint32_t makeForceContributor(apf::MeshEntity* e);
	virtual uint32_t makeConstraint(apf::MeshEntity* e);
	virtual uint32_t recover();

private:
	uint32_t integration_order;
	uint32_t polynomial_order;
	apf::Mesh* m;
	apf::Field* field;
	apf::Matrix< 3,3> D;
	apf::Numbering* nodeNums;
	apf::Numbering* faceNums;
	AlgebraicSystem* linsys;

};

/*expose this internal method for testing purposes*/
apf::Matrix< 3,3 > buildD(double E, double Nu, bool use_plane_stress);

#endif
