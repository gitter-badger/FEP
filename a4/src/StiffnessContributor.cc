#include <apfMesh.h>

#include "StiffnessContributor.h"

StiffnessContributor::StiffnessContributor(
	apf::Field* f,
	uint32_t integrate_order) : apf::Integrator(integrate_order), field(f)
{
  ndims = apf::getMesh(f)->getDimension();
}

void StiffnessContributor::inElement(apf::MeshElement* me)
{
	/*create the Field Element from the MeshElement*/
	this->field_element = apf::createElement(this->field, me);
	/*determine the size of the force matrices*/
	this->ndofs = this->ndims * apf::countNodes(this->field_element);
	this->fe.setSize(this->ndofs);
	this->ke.setSize(this->ndofs, this->ndofs);
	/*zero the internal matricies*/
	this->ke.zero();
	this->fe.zero();

}

void StiffnessContributor::outElement()
{
	/*perform clean up and destroy specific field element*/
	apf::destroyElement(this->field_element);
}

void StiffnessContributor::atPoint(apf::Vector3 const& p, double w, double dV)
{
	/**/
	std::cout << "p: " << p << std::endl;
	std::cout << "weight: " << w << " dV: " << dV << std::endl;

	apf::MeshElement* me = apf::getMeshElement(this->field_element);
	apf::Vector3 x;
	apf::mapLocalToGlobal(me, p, x);
	std::cout << "global x: " << x << std::endl;
}
