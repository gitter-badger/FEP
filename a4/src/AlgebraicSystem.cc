#include <stdexcept>
#include <climits>
#include <assert.h>

#include <petscmat.h>

#include "AlgebraicSystem.h"

AlgebraicSystem::AlgebraicSystem(std::size_t ngd)
	: nGlobalDOFs(ngd)
{
	this->_allow_assembly = false;
	/*preallocate the map since its keys will run from [0, nGlobalDOFs-1]*/
	for(std::size_t ii = 0; ii < this->nGlobalDOFs; ++ii) {
		/* -1 is valid intializer for unsigned types since this puts us at
		* UINT_MAX for what ever type we are using. Using this value as 
		* sentinal value for unitialized is valid, since this number
		* represents the indice on a square matrix. Since PETSc,
		* which will be handling this matrix, has never done more than
		* 500 million rows in a matrix before on super computers, it is
		* reasonable to assume that we will never be in danger of overflow
		* for this value for the next 20 years*/
		this->masks.insert(
			std::map< size_t, uint64_t >::value_type(ii, UNMAPPED_VALUE));
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
		throw std::invalid_argument(
			"cannot add constraint after assembly has started");
	}
	std::size_t input_size = fixed.size();
	/*check that there are at least as many mappings as there are elements
	* in fixed that will be mapped. Otherwise we will walk off end of 
	* buffer*/
	if(input_size > node_mapping.size()) {
		throw std::range_error(
			"local force vector rows size exceeded mapping size");
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
		/*add the proscribed displacements to our known_d(isplacement) store*/
		this->known_d.push_back(fixed[ii]);
	}

}

void AlgebraicSystem::beginAssembly() {
	_allow_assembly = true;

	/*now compute the remaining degrees of freedom*/
	std::size_t ndogs = this->known_d.size();
	std::size_t n_eq = this->nGlobalDOFs - ndogs;
	/*go throught the map and map the free degrees of freedom
	* to indicies that will be used with Petsc matricies and vectors
	* to solve the remaining degrees of freedom*/
	uint64_t n_eq_counter = 0;
	for(std::size_t ii = 0; ii < this->nGlobalDOFs; ++ii) {
		if(this->masks[ii] == UNMAPPED_VALUE) {
			/*clear the most significant bit to indicate that this is
			* a free dof */
			this->masks[ii] = n_eq_counter & (~(KNOWN_DOF_MASK));
			++n_eq_counter;
		} else {
			/*this dof was already labled*/
		}
	}
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

	ierr = MatCreateSBAIJ(PETSC_COMM_WORLD, 1, n_eq, n_eq, n_eq, n_eq,
		300, PETSC_NULL, 300, PETSC_NULL, &(this->K));
	ierr = MatSetOption(this->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetOption(this->K, MAT_SYMMETRIC, PETSC_TRUE);
	/*this part will ignore any insertions in the lower triangular, solving
	* all of my problems and making the assembly code 10x less complex */
	ierr = MatSetOption(this->K, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
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
	if(nLocalDOFs > this->nGlobalDOFs) {
		throw std::range_error("local mapping size exceeded global size");
	}
	/*this will accept a non square matrix as long as its mappings fit*/
	/*check that every row will have a mapping*/
	if((nrows > this->nGlobalDOFs) || (nrows > nLocalDOFs)) {
		throw std::range_error("local rows size exceeded mapping size");
	}
	/*check that every column will have mapping and will fit into the
	* global matrix*/
	if((ncol > this->nGlobalDOFs) || (ncol > nLocalDOFs)) {
		throw std::range_error("local column size exceeded mapping size");
	}
	/*it will be cheaper to allocate the full possible size for the
	* global indice vector than to compute exact size */
	PetscInt *idxM = new PetscInt[nrows];
	PetscInt *idxN = new PetscInt[ncol];
	PetscInt *ix = new PetscInt[nrows];

	/*determine size of the assembly matrix & any force vector contributions*/
	
	/*cache the local indicies we will need*/
	std::vector< std::size_t > keep_rows, keep_cols, add_to_f_col;
	keep_rows.clear();
	keep_cols.clear();
	add_to_f_col.clear();
	std::vector< double > local_displacements;
	local_displacements.clear();

	/*find the indices of the element stiffness we will keep*/
	PetscInt idxM_size = 0;

	for(std::size_t ii = 0; ii < nrows; ++ii){
		/*get the reduced global index using the map*/
		int tmp_global_index = node_mapping[ii];
		if(tmp_global_index < 0){
			throw std::logic_error(
				"cannot have negative indices for global dofs");
		}
		/*convert to unsigned type using promotion*/
		std::size_t tmp_key = tmp_global_index;
		uint64_t tmp_val = this->masks[tmp_key];
		/*test if dof is fixed or not*/
		if(!(tmp_val & KNOWN_DOF_MASK)){
			/*extract the 32 bit number*/
			idxM[idxM_size] = (PetscInt) (tmp_val & SAMPLE_LOW_32_BITS);
			//std::cout << "ii: " << ii << " maps to: " << idxM[idxM_size] << std::endl;
			++idxM_size;
			keep_rows.push_back(ii);
		} else {
			//std::cout << "dof already mapped" << std::endl;
		}
	}

	/*find the columns of the element stiffness matrix we will keep and 
	* the columns that will contribute to the force vector*/
	PetscInt idxN_size = 0;
	PetscInt ni = 0;
	for(std::size_t ii = 0; ii < ncol; ++ii) {
		/*get the reduced global index using the map*/
		int tmp_global_index = node_mapping[ii];
		if(tmp_global_index < 0){
			throw std::logic_error(
				"cannot have negative indices for global dofs");
		}
		/*convert to unsigned type using promotion*/
		std::size_t tmp_key = tmp_global_index;
		uint64_t tmp_val = this->masks[tmp_key];
		/*test if dof is fixed or not*/
		if(!(tmp_val & KNOWN_DOF_MASK)){
			/*high bit not set means it is a free degree*/
			/*extract the 32 bit number*/
			idxN[idxN_size] = (PetscInt) (tmp_val & SAMPLE_LOW_32_BITS);
			//std::cout << "ii: " << ii << " maps to: " << idxM[idxM_size] << std::endl;
			++idxN_size;
			keep_cols.push_back(ii);
		} else {
			/*high bit IS set, so it is fixed*/

			ix[ni] = (PetscInt) (tmp_val & SAMPLE_LOW_32_BITS);
			add_to_f_col.push_back(ii);
			local_displacements.push_back(this->known_d[(tmp_val & SAMPLE_LOW_32_BITS)]);
			
			std::cout << "col contributes to force vector: " << ii << " mapping: " << node_mapping[ii] <<  std::endl;
			
			++ni;
		}
	}
	assert(idxM_size == keep_rows.size());
	assert(idxN_size == keep_cols.size());
	assert(ni == local_displacements.size());
	assert(ni == add_to_f_col.size());

	/*allocate a two dimensional storage for stiffness matrix 
	* using PetscInts which by the default options chosen elsewhere
	* should be row-oriented*/
	PetscScalar *values = new PetscScalar[idxM_size * idxN_size];
	PetscScalar *y = new PetscScalar[idxM_size];

	/*loop over element stiffness matricies and copy contents into 
	* the petsc approved formats*/

	for(std::size_t ii = 0; ii < keep_rows.size(); ++ii){
		/*skip all rows that are already determined*/
		for(std::size_t jj = 0; jj < keep_cols.size(); ++jj){
			/*get the indicies we want in the local stiffness matrix*/
			uint64_t kk = keep_rows[ii];
			uint64_t ll = keep_cols[jj];
			/*load the petsc scalar values which are known to be not zero
			* by virtue of check above from the double*/
			PetscScalar tmp = ke(kk,ll);
			std::size_t computed_index = ii * idxM_size + jj;
			values[computed_index] = tmp;
		}
		/*move the stiffness contribution for the fixed dofs to the force
		* vector side*/
		for(std::size_t jj = 0; jj < add_to_f_col.size(); ++jj) {
			uint64_t kk = keep_rows[ii];
			uint64_t ll = add_to_f_col[jj];
			/*load the petsc scalar values which are known to be not zero
			* by virtue of check above from the double*/
			PetscScalar tmp = ke(kk,ll);
			tmp *= ( -1 * local_displacements[jj]);
			y[jj] = tmp;
		}
	}
	/*now that all of the above has been computed*/
	MatSetValues(this->K, idxM_size, idxM, idxN_size, idxN, values, ADD_VALUES);
	VecSetValues(this->F, ni, ix, y, ADD_VALUES);

	delete[] y;
	delete[] ix;
	delete[] idxM;
	delete[] idxN;
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
	if(nLocalDOFs > this->nGlobalDOFs) {
		throw std::range_error("local mapping size exceeded global size");
	}
	/*check that every row will have a mapping*/
	if((nrows > this->nGlobalDOFs) || (nrows > nLocalDOFs)) {
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
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	MatView(this->K, PETSC_VIEWER_STDOUT_WORLD);
	VecView(this->F, PETSC_VIEWER_STDOUT_WORLD);


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
