#ifndef ALGEBRAIC_SYSTEM
#define ALGEBRAIC_SYSTEM 1

#include <vector>
#include <map>
#include <stdint.h>

#include <apfDynamicMatrix.h>

#include <petscksp.h>

#include "BandedSymmetricMatrix.h"

class AlgebraicSystem
{
public:
	AlgebraicSystem(std::size_t global_n_dofs);
	~AlgebraicSystem();

	void assemble(
		apf::DynamicMatrix const& ke,
		apf::NewArray<int> const& node_mapping,
		uint32_t nLocalDOFs);
	void assemble(
		std::vector<double> const& fe,
		apf::NewArray<int> const& node_mapping,
		uint32_t nLocalDOFs);

private:
	Mat K;
    Vec d;
    Vec F;
    KSP solver;

	std::size_t ndofs;
	/*we keep track of which degrees of freedom are
	* fixed by boundary constraints first*/
	std::map<std::size_t, uint64_t> mask;
};

#endif
