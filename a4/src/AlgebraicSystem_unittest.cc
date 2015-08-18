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
	apf::DynamicMatrix ke_1;
	apf::NewArray<int> mapping_1;
	uint32_t ke_1_size;
	/*matrix designed to fit on lower diagonal of global matrix*/
	apf::DynamicMatrix ke_2;
	apf::NewArray<int> mapping_2;
	uint32_t ke_2_size;

	virtual void SetUp() {
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
		mapping_1[2] = 9;
		mapping_1[3] = 11;
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
		mapping_2[0] = 8;
		mapping_2[1] = 9;
		mapping_2[2] = 10;
		mapping_2[3] = 11;
	}

	virtual void TearDown() {

	}

};

TEST_F(AlgebraicSystemTest, AssembleInsertionTest) {
	/*use nGlobalDOFs, no boundary conditions for this one*/
	std::size_t nGlobalDOFs = 12;
	AlgebraicSystem linsys(nGlobalDOFs);

}

TEST_F(AlgebraicSystemTest, IgnoreLowerDiagonals) {
	/*use nGlobalDOFs, no boundary conditions for this one*/
	std::size_t nGlobalDOFs = 12;
	AlgebraicSystem linsys(nGlobalDOFs);

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
				EXPECT_EQ(ke_2(ii, kk), vals[jj]);
				
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
	std::size_t nGlobalDOFs = 12;
	AlgebraicSystem linsys(nGlobalDOFs);
	/*try to call assemble before calling beginAssembly
	* here we test for failure in form of execption*/
	EXPECT_THROW(linsys.assemble(ke_1, mapping_1, ke_1_size), std::invalid_argument);
}

TEST_F(AlgebraicSystemTest, AllowCorrectAssembly) {
	/*use nGlobalDOFs, no boundary conditions for this one*/
	std::size_t nGlobalDOFs = 12;
	AlgebraicSystem linsys(nGlobalDOFs);
	/*try to call assemble after calling beginAssembly
	* here we do not want an exception*/
	linsys.beginAssembly();
	EXPECT_NO_THROW(linsys.assemble(ke_1, mapping_1, ke_1_size));
}

	// PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

	// MatView(linsys.K, PETSC_VIEWER_STDOUT_WORLD);
	// VecView(linsys.F, PETSC_VIEWER_STDOUT_WORLD);
