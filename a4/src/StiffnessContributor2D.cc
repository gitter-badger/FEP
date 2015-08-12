#include <apfMesh.h>

#include "StiffnessContributor2D.h"

StiffnessContributor2D::StiffnessContributor2D(
	apf::Field* f,
	apf::Matrix< 3,3 > D_input,
	uint32_t integrate_order) : apf::Integrator(integrate_order), field(f), D(D_input)
{
  ndims = apf::getMesh(f)->getDimension();
}

void StiffnessContributor2D::inElement(apf::MeshElement* me)
{
	/*create the Field Element from the MeshElement*/
	this->field_element = apf::createElement(this->field, me);
	/*determine the size of the force matrices*/
	this->nnodes = apf::countNodes(this->field_element);
	this->ndofs = this->ndims * apf::countNodes(this->field_element);
	this->ke.setSize(this->ndofs, this->ndofs);
	/*zero the internal matricies*/
	this->ke.zero();
}

void StiffnessContributor2D::outElement()
{
	/*perform clean up and destroy specific field element*/
	apf::destroyElement(this->field_element);
}

void StiffnessContributor2D::atPoint(apf::Vector3 const& p, double w, double dV)
{
	apf::NewArray<double> shape_val;
  	apf::getShapeValues(this->field_element, p, shape_val);

 	apf::NewArray<apf::Vector3> gradShape;
 	apf::getShapeGrads(this->field_element, p, gradShape);

 	apf::Vector3 x;
  	apf::MeshElement* me = apf::getMeshElement(this->field_element);
  	apf::mapLocalToGlobal(me, p, x);

	//std::cout << this->D << std::endl;

	/*explicitly instandiate the size of this array*/
	/*this goes by columns, rows*/
	apf::NewArray< apf::Matrix< 3,2 > > B(this->nnodes);
	/*construct each of the nnodes shape function matricies*/
	uint32_t ii;
	for(ii = 0; ii < this->nnodes; ++ii) {
		B[ii][0][0] = gradShape[ii][0];
		B[ii][0][1] = 0.0;
		B[ii][1][0] = 0.0;
		B[ii][1][1] = gradShape[ii][1];
		B[ii][2][0] = gradShape[ii][1];
		B[ii][2][1] = gradShape[ii][0];
		std::cout << std::endl << B[ii] << std::endl;
		// returns the transpose of matrix; apf::transpose(B[ii]);
	}
	std::cout << "-------------------------------" << std::endl;

}
void StiffnessContributor2D::computeBtranspose(apf::Vector3 const& grad, apf::Matrix< 2,3 > & B_T)
{
	B_T[0][0] = grad[0];
}

void StiffnessContributor2D::computeB(apf::Vector3 const& grad, apf::Matrix< 3,2 > & B)
{

}

