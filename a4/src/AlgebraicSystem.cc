#include <stdexcept>
#include <climits>

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
		this->masks.insert(std::map< size_t, uint64_t >::value_type(ii, -1));
	}
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

void AlgebraicSystem::beginAssembly() {
	_allow_assembly = true;


	/*now compute the remaining degrees of freedom*/
	std::size_t ndogs = this->known_d.size();
	std::size_t n_eq = this->ndofs - ndogs;
	/*we construct the stiffness matrix only for degress of freedom
	* which have not been fixed, 
	* WARNING: this program is intended to run on a single process
	* there for the local size of all petsc variables is set to
	* the same as the global size.*/
	VecCreateMPI(PETSC_COMM_WORLD, n_eq, n_eq, &(this->F));
	std::cout << "here" << std::endl;
	VecSetOption(this->F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
		std::cout << "here" << std::endl;

	MatCreateAIJ(PETSC_COMM_WORLD, n_eq, n_eq, n_eq, n_eq,
		300, PETSC_NULL, 300, PETSC_NULL, &(this->K));
		std::cout << "here" << std::endl;

	MatSetOption(this->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		std::cout << "here" << std::endl;

	KSPCreate(PETSC_COMM_WORLD, &(this->solver));
		std::cout << "here" << std::endl;

	KSPSetTolerances(this->solver, 1.0e-8, 1.0e-8, PETSC_DEFAULT, 100);
		std::cout << "here" << std::endl;

	VecDuplicate(this->F, &(this->d));
		std::cout << "here" << std::endl;
}

void AlgebraicSystem::addBoundaryConstraint(
	std::vector<double> fixed,
	apf::NewArray<int> const& node_mapping,
	uint32_t size)
{
	/*do not add a constraint after assembly has begun*/
	if(this->_allow_assembly == true) {
		throw std::invalid_argument("cannot add constraint after assembly has started");
	}

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

	// for(std::size_t ii = 0; ii < nrows; ++ii){
	// 	for(std::size_t jj = 0; jj < ncol; ++jj){
	// 		/*lookup where we put this in the global array*/
	// 		std::size_t kk = node_mapping[ii];
	// 		std::size_t ll = node_mapping[jj];
	// 		this->K(kk,ll) += ke(ii,jj);
	// 	}
	// }
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
