#ifndef ALGEBRAIC_SYSTEM
#define ALGEBRAIC_SYSTEM 1

#include <vector>
#include <apfDynamicMatrix.h>

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
	std::size_t ndofs;
	BandedSymmetricMatrix K;
	std::vector< double > F;

};

#endif
