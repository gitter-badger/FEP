#ifndef ALGEBRAIC_SYSTEM
#define ALGEBRAIC_SYSTEM 1

#include <vector>
#include <map>
#include <stdint.h>

#include <apfDynamicMatrix.h>

#include <petscksp.h>

#define KNOWN_DOF_MASK ((uint64_t)1<<63)
#define SAMPLE_LOW_32_BITS ((uint32_t)(0xFFFFFFFF))
#define UNMAPPED_VALUE ((uint64_t)-1)

class AlgebraicSystem
{
public:
	AlgebraicSystem(std::size_t nGlobalDOFs);
	~AlgebraicSystem();

	void addBoundaryConstraint(
		std::vector<double> const& fixed,
		std::vector<uint32_t> const& node_mapping);

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
	
	Mat K;
	Vec F;

private:
	bool _allow_assembly;
	/*the key is the index in the global numbering of dofs,
	* most significant bit of the value is 1 if the global dof
	* is fixed, with the lower 32 bits being the only valid indices
	* because of interaction with PetscInt not being a reliable 
	* size across platforms*/
	std::map< std::size_t,uint64_t > masks;
	std::vector< double > known_d;

    Vec d;

    KSP solver;

	std::size_t nGlobalDOFs;
};

#endif
