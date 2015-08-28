#include <stdint.h>
#include <stdio.h>
#include <iomanip>
#include <cmath>

#include <gtest/gtest.h>

#include <apf.h>
#include <apfMesh2.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>

#include "StiffnessContributor2D.h"
#include "MeshBuilder.h"
#include "ElasticAnalysis2D.h"

class StiffnessTest : public testing::Test,
	public ::testing::WithParamInterface<int>
{
protected:
	apf::Mesh2* mesh;
	MeshBuilder* mesh_builder;
	apf::Field* field;
	apf::Matrix< 3,3 > D;
	double E;
	double Nu;
	int integration_order;

	void changeMeshFromIndex(int index) {
		if(NULL != this->mesh) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
		int X_ELMS = 2;
		int Y_ELMS = 1;
		switch(index) {
			case 0:
				mesh_builder->build2DRectQuadMesh(this->mesh, 2, 1, 0.0, 0.0, 2.0, 1.0);
				apf::changeMeshShape(this->mesh, apf::getLagrange(1));
				break;
			case 1:
				this->mesh_builder->build2DRectQuadMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(this->mesh, apf::getLagrange(2));
				break;
			case 2:
				this->mesh_builder->build2DRectTriMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				break;
			case 3:
				this->mesh_builder->build2DRectTriMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(this->mesh, apf::getLagrange(2));
				break;
			case 4:
				this->mesh_builder->build2DRectQuadMesh(this->mesh, X_ELMS, Y_ELMS,
					0.0, 0.0, 2.0, 2.0);
				apf::changeMeshShape(this->mesh, apf::getSerendipity());
			default:
				/*mesh should stay null*/
				break;
		}
		bool use_plane_stress = true;
		this->D = buildD(this->E, this->Nu, use_plane_stress);
		this->field = createField(this->mesh, "dummy", apf::VECTOR, this->mesh->getShape());
		apf::zeroField(this->field); /*must zero the field to force writing of data*/

	}

	virtual void SetUp() {
		mesh_builder = new MeshBuilder();
		mesh = NULL;
		this->E = 1e8;
		this->Nu = 0.35;
		this->integration_order = 4;
	}

	virtual void TearDown() {
		if(mesh != NULL) {
			mesh->destroyNative();
			apf::destroyMesh(mesh);
		}
		delete mesh_builder;	
	}

};

INSTANTIATE_TEST_CASE_P(DifferentMeshOrders, StiffnessTest,
	::testing::Values(0,1,2,3,4));

TEST_P(StiffnessTest, StiffnessIsSymmetric) {
	/*use a linear quad*/
	changeMeshFromIndex(GetParam());

	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = this->mesh->begin(2);
	while((e = this->mesh->iterate(it))) {
		/*compute each of the nodal submatrices by different method, and
		* do not explictly take advantage of symmetry*/
		StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

		apf::MeshElement* me = apf::createMeshElement(this->mesh, e);
		stiff.process(me);
		apf::destroyMeshElement(me);

		// std::cout << "======================" << std::endl;
		// std::cout << stiff.ke << std::endl;
		// std::cout << "======================" << std::endl;


		apf::Element* f_elm = apf::createElement(this->field, me);
		uint32_t nnodes = apf::countNodes(f_elm);
		apf::destroyElement(f_elm);

		uint32_t n_eqs = nnodes * this->mesh->getDimension();
		/*check that stiffness are symmetric*/
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = ii; jj < n_eqs; ++jj) {
				//EXPECT_FLOAT_EQ(stiff.ke(ii,jj), stiff.ke(jj,ii));
			}
		}
		/*visualization for the symmetry of the matrix*/
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = 0; jj < n_eqs; ++jj) {
				if(fabs(stiff.ke(ii,jj) - stiff.ke(jj,ii)) > 1e-13) {
					std::cout << "X ";
				} else {
					std::cout << "0 ";
				}
			}
			std::cout << std::endl;
		}
	}
	this->mesh->end(it);

}

TEST_P(StiffnessTest, CheckStiffnessMatrix) {
	/*use a linear quad*/
	changeMeshFromIndex(GetParam());

	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = this->mesh->begin(2);
	while((e = this->mesh->iterate(it))) {
		/*compute each of the nodal submatrices by different method, and
		* do not explictly take advantage of symmetry*/
		StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

		apf::MeshElement* me = apf::createMeshElement(this->mesh, e);

		stiff.process(me);

		apf::Element* f_elm = apf::createElement(this->field, me);

		uint32_t nnodes = apf::countNodes(f_elm);
		uint32_t n_eqs = nnodes * this->mesh->getDimension();

		apf::DynamicMatrix alt_ke;

		alt_ke.setSize(n_eqs,n_eqs);
		alt_ke.zero();

		int entity_type = this->mesh->getType(e);
		uint32_t nIntPoints = apf::countGaussPoints(entity_type, this->integration_order);

		for(uint32_t ii = 0; ii < nIntPoints; ++ii) {
			apf::Vector3 p(0.0,0.0,0.0);
			apf::getGaussPoint(entity_type, this->integration_order, ii, p);
			double dV = me->getDV(p);

			double weight = apf::getIntWeight(me, this->integration_order, ii);

			apf::NewArray<apf::Vector3> gradShape;

			apf::getShapeGrads(f_elm, p, gradShape);

			apf::Matrix<2,2> nodal_stiffness;

			for(uint32_t A = 0; A < nnodes; ++A) {
				for(uint32_t B = 0; B < nnodes; ++B) {
					// if(A != B) {
					// 	continue;
					// }
					/* nodal_stiffness_11
					* in reduced form we have N_A,x * D_11 * N_B,x
											+ N_A,y * D_33 * N_B,y*/
					double tmp = gradShape[A][0] * this->D[0][0] * gradShape[B][0];
					tmp += gradShape[A][1] * this->D[2][2] * gradShape[B][1];
					tmp *= (dV * weight);
			//		alt_ke(2*A,2*B) += tmp;
					
					/* nodal_stiffness_12
					* we have N_A,x * D_12 * N_B,y
							+ N_A,y * D_33 * N_B,x*/
					tmp = 0.0;
					tmp = gradShape[A][0] * this->D[0][1] * gradShape[B][1];
					tmp += gradShape[A][1] * this->D[2][2] * gradShape[B][0];
					tmp *= (dV * weight);
			//		alt_ke(2*A, 2*B + 1) += tmp;

					/* nodal_stiffness_21
					*we have N_A,y * D_12 * N_B,x
							+N_A,x * D_33 * N_B,y*/
					tmp = gradShape[A][1] * this->D[0][1] * gradShape[B][0];
					tmp += gradShape[A][0] * this->D[2][2] * gradShape[B][1];
					tmp *= (dV * weight);
			//		alt_ke(2*A + 1, 2*B) += tmp;

					/* nodal_stiffness_22
					*we have N_A,y * D_22 * N_B,y
							+N_A,x * D_33 * N_B,x*/
					tmp = gradShape[A][1] * this->D[1][1] * gradShape[B][1];
					tmp += gradShape[A][0] * this->D[2][2] * gradShape[B][0];
					tmp *= (dV * weight);
					alt_ke(2*A + 1, 2*B + 1) += tmp;

					//alt_ke(2*A, 2*B) += 1e10;
				//	alt_ke(2*A, 2*B + 1) += 1e10;
					//alt_ke(2*A + 1, 2*B) += 1e10;
					//alt_ke(2*A + 1, 2*B + 1) += 1e10;
				}
			}
		}
		/*check that altenate calculation of stiffness is symmetric*/
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = ii+1; jj < n_eqs; ++jj) {
				//EXPECT_FLOAT_EQ(alt_ke(ii,jj), alt_ke(jj,ii));
			}
		}
		std::cout << "=======================================" << std::endl;
		std::cout << "=			Reference Matrix 			=" << std::endl;
		std::cout << alt_ke << std::endl;
		std::cout << "=======================================" << std::endl << std::endl;

		/*visualization for the symmetry of the matrix*/
		std::cout << std::setprecision(16);
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = 0; jj < n_eqs; ++jj) {
				if(fabs(alt_ke(ii,jj) - alt_ke(jj,ii)) > fabs(1e-10 * alt_ke(ii,jj))) {
					// std::cout << "X ";
					/*absolute error*/
					if(fabs(alt_ke(ii,jj) - alt_ke(jj,ii)) > 5e-9) {
						std::cout << fabs(alt_ke(ii,jj) - alt_ke(jj,ii)) << " ";
					} else {
						std::cout << "0 ";
					}

				} else {
					std::cout << "0 ";
				}
			}
			std::cout << std::endl;
		}
		std::cout << "-----------------------------------" << std::endl;
		// std::cout << "=======================================" << std::endl;
		// std::cout << "=			Difference between two		=" << std::endl;
		// std::cout << "=======================================" << std::endl;
		// for(uint32_t ii = 0; ii < n_eqs; ++ii) {
		// 	for(uint32_t jj = 0; jj < n_eqs; ++jj) {
		// 		if(fabs(alt_ke(ii,jj) - stiff.ke(ii,jj)) > 1e-13) {
		// 			std::cout << "X ";
		// 		} else {
		// 			std::cout << "0 ";
		// 		}
		// 	}
		// 	std::cout << std::endl;
		// }



		/*now check that we computed the same stiffness matrix*/
		for(uint32_t ii = 0; ii < n_eqs; ++ii) {
			for(uint32_t jj = 0; jj < n_eqs; ++jj) {
				//EXPECT_FLOAT_EQ(alt_ke(ii,jj), stiff.ke(ii,jj));
			}
		}
		apf::destroyElement(f_elm);
		apf::destroyMeshElement(me);
	}
	this->mesh->end(it);
}

class MyUserFunction : public apf::Function
{
public:
	MyUserFunction(apf::Field* f, apf::Mesh* m) : node_field(f), mesh(m) {}

	void eval(apf::MeshEntity* e, double *result) {
		if(NULL != node_field) {
			apf::MeshElement* me = apf::createMeshElement(this->mesh, e);
			apf::Element* f_elm = apf::createElement(this->node_field, me);
			/*find the location of the single node on entity*/
			apf::Vector3 tmp_vec;
			apf::getVector(this->node_field, e, 0, tmp_vec);

			apf::Vector3 x_vec;
			apf::mapLocalToGlobal(me, tmp_vec, x_vec);
			*result = (2*x_vec[0] + 5*x_vec[1]);
		} else {
			std::cout << "failed to set field" << std::endl;
			*result = -1.0;
		}

	}

protected:
	apf::Field* node_field;
	apf::Mesh* mesh;

};

TEST_F(StiffnessTest, ScratchPad) {
	apf::Mesh2* m = NULL;
	this->mesh_builder->buildBatmanElementMesh(m);

	apf::MeshIterator* it;
	apf::MeshEntity* e;
	it = m->begin(2);
	e = m->iterate(it);
	m->end(it);

	apf::Field* node_f = apf::createField(m, "nodeField", apf::SCALAR, m->getShape());
	apf::zeroField(node_f);
	MyUserFunction test_fnc(node_f, m);
	apf::Field* test_f = apf::createUserField(m, "testingField", apf::SCALAR, m->getShape(), &test_fnc);

	apf::MeshElement* me = apf::createMeshElement(m, e);
	apf::Element* f_elm = apf::createElement(test_f, me);

	/*pick a random spot inside the element coordinates*/
	apf::Vector3 sample_points[8] =
	{apf::Vector3(-0.1241,    0.4889, 0.0),
	 apf::Vector3( 0.4897,    0.0347, 0.0),
	 apf::Vector3( 0.4090,    0.7269, 0.0),
	 apf::Vector3( 0.4172,   -0.3034, 0.0),
	 apf::Vector3( 0.6715,    0.2939, 0.0),
	 apf::Vector3(-0.2075,   -0.7873, 0.0),
	 apf::Vector3( 0.7172,    0.8884, 0.0),
	 apf::Vector3( 0.6302,   -0.1471, 0.0)};
	for(uint32_t ii = 0; ii < 8; ++ii) {
		apf::Vector3 tmp_grad(0.0, 0.0, 0.0);
		apf::getGrad(f_elm, sample_points[ii], tmp_grad);
		EXPECT_FLOAT_EQ(2.0, tmp_grad[0]);
		EXPECT_FLOAT_EQ(5.0, tmp_grad[1]);
	} 
	apf::writeVtkFiles("batman_elm", m);
}
