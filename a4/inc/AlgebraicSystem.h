#ifndef ALGEBRAIC_SYSTEM
#define ALGEBRAIC_SYSTEM 1

#include <vector>
#include <map>
#include <stdint.h>
#include <apfDynamicMatrix.h>


#include "BandedMaskedMatrix.h"

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

	BandedMaskedMatrix K;
	std::vector< double > F;
private:
	std::size_t ndofs;
	/*we keep track of which degrees of freedom are
	* fixed by boundary constraints first*/
	std::map<std::size_t, uint64_t> mask;


};

#endif
