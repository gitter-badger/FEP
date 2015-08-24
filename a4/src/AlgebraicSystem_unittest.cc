#include <stdint.h>
#include <stdio.h>
#include <iomanip>

#include <apf.h>
#include <apfDynamicMatrix.h>

#include <gtest/gtest.h>

#include <petscmat.h>

#include "AlgebraicSystem.h"

class AlgebraicSystemTest : public testing::Test
{
protected:
	std::size_t nGlobalDOFs;


	apf::DynamicMatrix ke_1;
	apf::NewArray<int> mapping_1;
	uint32_t ke_1_size;
	/*matrix designed to fit on lower diagonal of global matrix*/
	apf::DynamicMatrix ke_2;
	apf::NewArray<int> mapping_2;
	uint32_t ke_2_size;
	/*fixture for some displacements testing*/
	std::vector<double> local_disp;
	std::vector<uint64_t> local_mapping;


	virtual void SetUp() {
		nGlobalDOFs = 12; /*must pick avalue above 4*/
		/*================================*/
		ke_1_size = 4;
		mapping_1.allocate(ke_1_size);	
		ke_1.setSize(ke_1_size, ke_1_size);
		for(std::size_t ii = 0; ii < ke_1_size; ++ii) {
			for(std::size_t jj = 0; jj < ke_1_size; ++jj) {
				if(ii == jj) {
					ke_1(ii, jj) = ii * jj;
				} else {
					ke_1(ii, jj) = 0.0;
				}
			}
		}
		/*this structure will spread out the submatrix
		* all over the larger matrix*/
		mapping_1[0] = 1;
		mapping_1[1] = 3;
		mapping_1[2] = nGlobalDOFs-3;
		mapping_1[3] = nGlobalDOFs-1;
		/*==================================*/
		ke_2_size = 4;
		mapping_2.allocate(ke_2_size);	
		ke_2.setSize(ke_2_size, ke_2_size);
		for(std::size_t ii = 0; ii < ke_2_size; ++ii) {
			for(std::size_t jj = 0; jj < ke_2_size; ++jj) {
				if(ii == jj) {
					/*on the main diagonal*/
					ke_2(ii, jj) = 10.0 * ii + jj + 1;
				} else if(ii > jj) {
					ke_2(ii,jj) = 7.0;
				} else {
					ke_2(ii, jj) = 0.0;
				}
			}
		}
		/*this structure will spread out the submatrix
		* all over the larger matrix*/
		mapping_2[0] = nGlobalDOFs-4;
		mapping_2[1] = nGlobalDOFs-3;
		mapping_2[2] = nGlobalDOFs-2;
		mapping_2[3] = nGlobalDOFs-1;
		/*================================*/
		local_disp.clear();
		local_disp.push_back(6.0);
		local_disp.push_back(8.0);
		local_disp.push_back(7.0);
		local_mapping.clear();
		local_mapping.push_back(nGlobalDOFs-4);
		local_mapping.push_back(nGlobalDOFs-1);
		local_mapping.push_back(0);
	}

	virtual void TearDown() {

	}

};

TEST_F(AlgebraicSystemTest, AssembleInsertionTest) {
	/*in this test we check that insertions of stiffness matrix
	* are additive, that is if we insert the same local stiffness 
	* matrix into the global matrix that the contributions are
	* additive. */

	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(nGlobalDOFs);

	linsys.beginAssembly();

	uint32_t nReps = 4;
	for(std::size_t ii = 0; ii < nReps; ++ii) {
		linsys.assemble(ke_1, mapping_1, ke_1_size);
	}
	linsys.synchronize();

	PetscInt ncols;
	const PetscInt *cols;
	const PetscScalar *vals;

	PetscInt row;
	/*check that the rows above the diagonal are correct*/
	for(std::size_t ii = 0; ii < ke_1_size; ++ii){
		row = mapping_1[ii];
		/*we only expect to see the upper triangular portion*/
		MatGetRowUpperTriangular(linsys.K);
		MatGetRow(linsys.K, row, &ncols, &cols, &vals);
		/*cross check with the initial submatrix*/
		for(std::size_t jj = 0; jj < ncols; ++jj) {
			std::size_t kk = -1;
			/*find the mapping for this column back to local
			* this is done by looping over every mapping until we find it
			* this is hacky but ok since the size is only 4*/
			for(std::size_t ll = 0; ll < ke_1_size; ++ll) {
				if(mapping_1[ll] == cols[jj]) {
					kk = ll;
					break; /*stop once we have found a mapping*/
				}
			}
			if(kk != -1) {
				EXPECT_TRUE(kk < ke_1_size);
				EXPECT_TRUE(ii < ke_1_size);
				/*the expected value is nReps times the local matrix*/
				EXPECT_FLOAT_EQ((ke_1(ii, kk) * nReps), vals[jj]);
				
			} else {
				/*we should always find mapping*/
				FAIL();
			}
		}
		MatRestoreRow(linsys.K, row, &ncols, &cols, &vals);
		MatRestoreRowUpperTriangular(linsys.K);
	}

}

TEST_F(AlgebraicSystemTest, IgnoreLowerDiagonals) {
	/*In this test we insert a local stiffness matrix
	* that has components above and below the main diagonal
	* of the global matrix, we ensure that these lower
	* contributions are not added to the upper triangular
	* portion of the global matrix by the underlying 
	* implementation. Basically we are checking that the
	* library provided matrix behaves the way we think
	* we are using it. */


	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(this->nGlobalDOFs);

	linsys.beginAssembly();
	linsys.assemble(ke_2, mapping_2, ke_2_size);
	linsys.synchronize();

	PetscInt ncols;
	const PetscInt *cols;
	const PetscScalar *vals;

	PetscInt row;
	/*check that the rows above the diagonal are correct*/
	for(std::size_t ii = 0; ii < ke_2_size; ++ii){
		row = mapping_2[ii];
		/*we only expect to see the upper triangular portion*/
		MatGetRowUpperTriangular(linsys.K);
		MatGetRow(linsys.K, row, &ncols, &cols, &vals);
		/*cross check with the initial submatrix*/
		for(std::size_t jj = 0; jj < ncols; ++jj) {
			std::size_t kk = -1;
			/*find the mapping for this column back to local
			* this is done by looping over every mapping until we find it
			* this is hacky but ok since the size is only 4*/
			for(std::size_t ll = 0; ll < ke_2_size; ++ll) {
				if(mapping_2[ll] == cols[jj]) {
					kk = ll;
					break; /*stop once we have found a mapping*/
				}
			}
			if(kk != -1) {
				EXPECT_TRUE(kk < ke_2_size);
				EXPECT_TRUE(ii < ke_2_size);
				/*only test when */
				EXPECT_FLOAT_EQ(ke_2(ii, kk), vals[jj]);
				
			} else {
				/*we should always find mapping*/
				FAIL();
			}
		}
		MatRestoreRow(linsys.K, row, &ncols, &cols, &vals);
		MatRestoreRowUpperTriangular(linsys.K);
	}
}

TEST_F(AlgebraicSystemTest, PreventPrematureAssembly) {
	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*try to call assemble before calling beginAssembly
	* here we test for failure in form of execption*/
	EXPECT_THROW(linsys.assemble(ke_1, mapping_1, ke_1_size), std::invalid_argument);
}

TEST_F(AlgebraicSystemTest, AllowCorrectAssembly) {
	/*use nGlobalDOFs, no boundary conditions for this one*/
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*try to call assemble after calling beginAssembly
	* here we do not want an exception*/
	linsys.beginAssembly();
	EXPECT_NO_THROW(linsys.assemble(ke_1, mapping_1, ke_1_size));
}

TEST_F(AlgebraicSystemTest, PreventMutltipleAssembly) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	EXPECT_NO_THROW(linsys.beginAssembly());
	EXPECT_THROW(linsys.beginAssembly(), std::logic_error);
}

TEST_F(AlgebraicSystemTest, RepeatedConstraints) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	EXPECT_NO_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping));
	/*over constraining the system is allowed, even if we add the same entity
	* twice, it is expected to overwrite the exisiting displacements. So to 
	* verify we change one of the displacements*/
	this->local_disp[0] += 27.0;
	EXPECT_NO_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping));
	/*the known displacements are private, so we must go ahead and solve
	* the empty system to then extract the solution vector in order to check
	* if the overwriting was sucessful. We are only solving a tiny system,
	* so this behaviour is acceptable*/
	linsys.beginAssembly();
	linsys.synchronize();
	linsys.solve();

	std::vector<double> solution;
	solution.clear();
	linsys.extractDisplacement(solution);
	EXPECT_EQ(nGlobalDOFs, solution.size());
	for(std::size_t ii = 0; ii < nGlobalDOFs; ++ii) {
		if(ii == this->local_mapping[0]) {
			EXPECT_EQ(this->local_disp[0], solution[ii]);
		} else if(ii == this->local_mapping[1]) {
			EXPECT_EQ(this->local_disp[1], solution[ii]);
		} else if(ii == this->local_mapping[2]) {
			EXPECT_EQ(this->local_disp[2], solution[ii]);
		} else {
			EXPECT_FLOAT_EQ(0.0, solution[ii]);
		}
	}
}

TEST_F(AlgebraicSystemTest, OutOfRangeContraint) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*corrupt by making this one clearly out of range*/
	this->local_mapping[1] = nGlobalDOFs+1;

	EXPECT_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping), std::invalid_argument);
}

TEST_F(AlgebraicSystemTest, ExtraMappingsAreIgnored) {
	AlgebraicSystem linsys(this->nGlobalDOFs);
	/*this one is clearly out of range but it will not be used since there are 
	* fewer displacements than there are mappings, so this last one will not be used*/
	local_mapping.push_back(nGlobalDOFs+1);
	/*quick test for consitency, this should always pass, only included
	* to prevent others from dicking with the test fixture and messing
	* this up*/
	ASSERT_EQ(this->local_mapping.size(), this->local_disp.size()+1);

	EXPECT_NO_THROW(linsys.addBoundaryConstraint(this->local_disp, this->local_mapping));
}

TEST_F(AlgebraicSystemTest, ExtractSolutionWithConstraints) {
	AlgebraicSystem linsys(this->nGlobalDOFs);

	linsys.addBoundaryConstraint(this->local_disp, this->local_mapping);
	linsys.beginAssembly();
	linsys.synchronize();
	linsys.solve();

	std::vector<double> solution;
	solution.clear();
	linsys.extractDisplacement(solution);
	EXPECT_EQ(nGlobalDOFs, solution.size());
	for(std::size_t ii = 0; ii < nGlobalDOFs; ++ii) {
		if(ii == this->local_mapping[0]) {
			EXPECT_EQ(this->local_disp[0], solution[ii]);
		} else if(ii == this->local_mapping[1]) {
			EXPECT_EQ(this->local_disp[1], solution[ii]);
		} else if(ii == this->local_mapping[2]) {
			EXPECT_EQ(this->local_disp[2], solution[ii]);
		} else {
			EXPECT_FLOAT_EQ(0.0, solution[ii]);
		}
	}
}

	// PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	// MatView(linsys.K, PETSC_VIEWER_STDOUT_WORLD);
	// VecView(linsys.F, PETSC_VIEWER_STDOUT_WORLD);
