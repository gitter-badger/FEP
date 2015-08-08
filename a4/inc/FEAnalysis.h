#ifndef FEANALYSIS_BASE
#define FEANALYSIS_BASE 1

#include <stdint.h>

#include <apf.h>
#include <apfMesh.h>

/*Abstract base class for all analysis codes*/
class FEAnalysis
{
public:
	virtual uint32_t setup() = 0;
	virtual uint32_t solve() = 0;
	virtual uint32_t makeStiffnessContributor(apf::MeshEntity* e) = 0;
	virtual uint32_t makeForceContributor(apf::MeshEntity* e) = 0;
	virtual uint32_t makeConstraint(apf::MeshEntity* e) = 0;
	virtual uint32_t recover() = 0;
};

#endif
