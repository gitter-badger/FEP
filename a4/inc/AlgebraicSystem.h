#ifndef ALGEBRAIC_SYSTEM
#define ALGEBRAIC_SYSTEM 1

#include <vector>
#include <map>
#include <stdint.h>

#include <apfDynamicMatrix.h>

#include <petscksp.h>

#include "BandedSymmetricMatrix.h"

#define KNOWN_DOF_MASK ((uint64_t)1<<63)

class AlgebraicSystem
{
public:
	AlgebraicSystem(std::size_t global_n_dofs);
	~AlgebraicSystem();

	void addBoundaryConstraint(
		std::vector<double> fixed,
		apf::NewArray<int> const& node_mapping,
		uint32_t size);

	void beginAssembly();

	void assemble(
		apf::DynamicMatrix const& ke,
		apf::NewArray<int> const& node_mapping,
		uint32_t size);
	void assemble(
		std::vector<double> const& fe,
		apf::NewArray<int> const& node_mapping,
		uint32_t size);

	PetscErrorCode synchronize();
	PetscErrorCode solve();

private:
	bool _allow_assembly;
	/*the key is the index in the global numbering of dofs,
	( most significant bit of the value is 1 if the global dof
	* is fixed, with the lower  */
	std::map< std::size_t,uint64_t > masks;
	std::vector< double > known_d;

	Mat K;
    Vec d;
    Vec F;
    KSP solver;

	std::size_t ndofs;
};

#endif
