#include <stdexcept>

#include "AlgebraicSystem.h"

AlgebraicSystem::AlgebraicSystem(std::size_t global_n_dofs) : ndofs(global_n_dofs)
{
	this->K.setSize(global_n_dofs);
	this->K.zero();
	this->F.clear();
	this->F.resize(global_n_dofs, 0.0);
}

AlgebraicSystem::~AlgebraicSystem()
{
	
}


void AlgebraicSystem::assemble(
		apf::DynamicMatrix const& ke,
		apf::NewArray<int> const& node_mapping,
		uint32_t nLocalDOFs)
{
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

	for(std::size_t ii = 0; ii < nrows; ++ii){
		for(std::size_t jj = 0; jj < ncol; ++jj){
			/*lookup where we put this in the global array*/
			std::size_t kk = node_mapping[ii];
			std::size_t ll = node_mapping[jj];
			this->K(kk,ll) += ke(ii,jj);
		}
	}
}

void AlgebraicSystem::assemble(
	std::vector<double> const& fe,
	apf::NewArray<int> const& node_mapping,
	uint32_t nLocalDOFs)
{
	
}
