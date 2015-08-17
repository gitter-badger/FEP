#include <stdexcept>
#include <climits>
#include <assert.h>

#include "AlgebraicSystem.h"

AlgebraicSystem::AlgebraicSystem(std::size_t global_n_dofs) : ndofs(global_n_dofs)
{
	this->_allow_assembly = false;
	/*preallocate the map since its keys will run from [0, ndofs-1]*/
	for(std::size_t ii = 0; ii < this->ndofs; ++ii) {
		/* -1 is valid intializer for unsigned types since this puts us at
		* UINT_MAX for what ever type we are using. Using this value as 
		* sentinal value for unitialized is valid, since this number
		* represents the indice on a square matrix. Since PETSc,
		* which will be handling this matrix, has never done more than
		* 500 million rows in a matrix before on super computers, it is
		* reasonable to assume that we will never be in danger of overflow
		* for this value for the next 20 years*/
		this->masks.insert(std::map< size_t, uint64_t >::value_type(ii, UNMAPPED_VALUE));
	}
	/*make sure there are zero elements in the known_d vector, following
	* code depends on it starting out as zero*/
	this->known_d.clear();
	assert(0 == this->known_d.size());
}

AlgebraicSystem::~AlgebraicSystem()
{
	/*only destroy the matricies if they were allocated*/
	if(true == this->_allow_assembly) {
		MatDestroy(&(this->K));
		VecDestroy(&(this->d));
		VecDestroy(&(this->F));
		KSPDestroy(&(this->solver));
	}
}

void AlgebraicSystem::addBoundaryConstraint(
	std::vector<double> const& fixed,
	std::vector<uint32_t> const& node_mapping)
{
	/*do not add a constraint after assembly has begun*/
	if(this->_allow_assembly == true) {
		throw std::invalid_argument("cannot add constraint after assembly has started");
	}
	std::size_t input_size = fixed.size();
	/*check that there are at least as many mappings as there are elements
	* in fixed that will be mapped. Otherwise we will walk off end of 
	* buffer*/
	if(input_size > node_mapping.size()) {
		throw std::range_error("local force vector rows size exceeded mapping size");
	}

	/*look up current index number of ndogs, we will start
	* inserting each fixed dof at this index growing the vector
	* of known displacements from here*/
	uint64_t cur_ndogs = this->known_d.size();

	for(std::size_t ii = 0; ii < input_size; ++ii) {
		/*convert values to fixed width, discarding extra precision
		* if so necessary by using Cs type promotion rules*/
		uint64_t key = node_mapping[ii];
		/*look up key an convert it */
		this->masks[key] = cur_ndogs | (KNOWN_DOF_MASK);
		printf("key %lx \n", this->masks[key]);
		++cur_ndogs;
		/*add the proscribed displacements to our known_displacement store*/
		this->known_d.push_back(fixed[ii]);
	}

}

void AlgebraicSystem::beginAssembly() {
	_allow_assembly = true;

	/*now compute the remaining degrees of freedom*/
	std::size_t ndogs = this->known_d.size();
	std::size_t n_eq = this->ndofs - ndogs;
	std::cout << "ndogs: " << ndogs << " n_eq: " << n_eq << std::endl;
	/*go throught the map and map the free degrees of freedom
	* to indicies that will be used with Petsc matricies and vectors
	* to solve the remaining degrees of freedom*/
	uint64_t n_eq_counter = 0;
	for(std::size_t ii = 0; ii < this->ndofs; ++ii) {
		if(this->masks[ii] == UNMAPPED_VALUE) {
			/*clear the most significant bit to indicate that this is
			* a free dof */
			this->masks[ii] = n_eq_counter & (~(KNOWN_DOF_MASK));
			++n_eq_counter;
		} else {
			/*this dof was already labled*/
			std::cout << "already fixed dof: " << ii << std::endl;
		}
	}
	std::cout << n_eq_counter << ": " << n_eq <<std::endl;
	/*we should have seen exactly n_eq free dofs*/
	assert(n_eq == n_eq_counter);

	/*we construct the stiffness matrix only for degress of freedom
	* which have not been fixed, 
	* WARNING: this program is intended to run on a single process
	* there for the local size of all petsc variables is set to
	* the same as the global size.*/
	PetscErrorCode ierr;
	ierr = VecCreateMPI(PETSC_COMM_WORLD, n_eq, n_eq, &(this->F));
	ierr = VecSetOption(this->F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
	ierr = MatCreateAIJ(PETSC_COMM_WORLD, n_eq, n_eq, n_eq, n_eq,
		300, PETSC_NULL, 300, PETSC_NULL, &(this->K));
	ierr = MatSetOption(this->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = KSPCreate(PETSC_COMM_WORLD, &(this->solver));
	ierr = KSPSetTolerances(this->solver, 1.0e-8, 1.0e-8, PETSC_DEFAULT, 100);
	ierr = VecDuplicate(this->F, &(this->d));
}

void AlgebraicSystem::assemble(
	apf::DynamicMatrix const& ke,
	apf::NewArray<int> const& node_mapping,
	uint32_t nLocalDOFs)
{
	/*CHECK if assembly is possible*/
	if(this->_allow_assembly == false) {
		throw std::invalid_argument("must call beginAssembly before assembly");
	}

	/*iterate over the input matrix and add each value
	* into the global matrix according to the mapping in */
	std::size_t nrows = ke.getRows();
	std::size_t ncol = ke.getColumns();
	/*sanity checks that we actually have a mapglobal_n_dofsping*/
	if(nLocalDOFs > this->ndofs) {
		throw std::range_error("local mapping size exceeded global size");
	}
	/*this will accept a non square matrix as long as its mappings fit*/
	/*check that every row will have a mapping*/
	if((nrows > this->ndofs) || (nrows > nLocalDOFs)) {
		throw std::range_error("local rows size exceeded mapping size");
	}
	/*check that every column will have mapping and will fit into the
	* global matrix*/
	if((ncol > this->ndofs) || (ncol > nLocalDOFs)) {
		throw std::range_error("local column size exceeded mapping size");
	}
	/*it will be cheaper to allocate the full possible size for the
	* global indice vector than to compute exact size */
	PetscInt *idxm = new PetscInt[nrows];
	PetscInt *idxn = new PetscInt[ncol];

	/*determine size of the assembly matrix and any force vector contributions*/

	std::size_t idxm_size = 0;
	std::size_t idxn_size = 0;

	/*allocate a two dimensional storage for stiffness matrix 
	* using PetscInts which by the default options chosen elsewhere
	* should be row-oriented*/
	PetscInt *values = new PetscInt[idxm_size * idxn_size] ;

	// for(std::size_t ii = 0; ii < nrows; ++ii){
	// 	for(std::size_t jj = 0; jj < ncol; ++jj){
	// 		/*lookup where we put this in the global array*/
	// 		std::size_t kk = node_mapping[ii];
	// 		std::size_t ll = node_mapping[jj];
	// 		this->K(kk,ll) += ke(ii,jj);
	// 	}
	// }

	delete[] idxm;
	delete[] idxn;
	delete[] values;
}

void AlgebraicSystem::assemble(
	std::vector<double> const& fe,
	apf::NewArray<int> const& node_mapping,
	uint32_t nLocalDOFs)
{
	/*iterate over the input vector and add each value
	* into the global force vector according to the mapping in */
	std::size_t nrows = fe.size();
	/*sanity checks that we actually have a mapglobal_n_dofsping*/
	if(nLocalDOFs > this->ndofs) {
		throw std::range_error("local mapping size exceeded global size");
	}
	/*check that every row will have a mapping*/
	if((nrows > this->ndofs) || (nrows > nLocalDOFs)) {
		throw std::range_error("local rows size exceeded mapping size");
	}


//	for(std::size_t ii = 0; ii < nrows; ++ii){
//		/*lookup where we put this in the global array*/
//		std::size_t kk = node_mapping[ii];
//		this->F[kk] += fe[ii];
//	}	
}

PetscErrorCode AlgebraicSystem::synchronize()
{
	PetscErrorCode ierr;
	ierr = VecAssemblyBegin(this->F);
  	CHKERRQ(ierr);
  	ierr = VecAssemblyEnd(this->F);
  	CHKERRQ(ierr);
  	ierr = MatAssemblyBegin(this->K, MAT_FINAL_ASSEMBLY);
  	CHKERRQ(ierr);
  	ierr = MatAssemblyEnd(this->K, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
  	return (PetscErrorCode)0;
}

PetscErrorCode AlgebraicSystem::solve()
{
	PetscErrorCode ierr;
	/*we use same matrix as preconditioner*/
	ierr = KSPSetOperators(this->solver, this->K, this->K);
	CHKERRQ(ierr);
  	ierr = KSPSetFromOptions(this->solver);
  	CHKERRQ(ierr);
  	ierr = KSPSolve(this->solver, this->F, this->d);
  	CHKERRQ(ierr);
  	/*zero is no error*/
  	return (PetscErrorCode)0;
}
