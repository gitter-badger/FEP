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
	std::cout << this->D << std::endl;
}
