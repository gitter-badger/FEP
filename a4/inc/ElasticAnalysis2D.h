#ifndef ELASTIC_ANALYSIS
#define ELASTIC_ANALYSIS 1

#include <apf.h>
#include <apfMesh.h>
#include <apfMatrix.h>

#include "FEAnalysis.h"

class ElasticAnalysis2D : FEAnalysis
{
public:
	ElasticAnalysis2D(apf::Mesh* m, uint32_t integration_order, double E, double nu);
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

};

#endif
