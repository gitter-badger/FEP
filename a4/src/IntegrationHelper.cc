#include <apfMesh.h>

#include "IntegrationHelper.h"

IntegrationHelper::IntegrationHelper(apf::Field* f, uint32_t order) :
  apf::Integrator(order),
  field(f)
{
  n_dimensions = apf::getMesh(f)->getDimension();
}

void IntegrationHelper::inElement(apf::MeshElement* me)
{
	/*create the Field Element from the MeshElement*/
	this->field_element = apf::createElement(this->field, me);
	/*determine the size of the force matrices*/
	this->n_dofs = apf::countNodes(this->field_element);
	this->f_e.setSize(this->n_dofs);
	this->k_e.setSize(this->n_dofs, this->n_dofs);
	/*zero the internal matricies*/
	uint32_t ii, jj;
	for(ii = 0; ii < this->n_dofs; ++ii) {
		this->f_e(ii) = 0.0;
		for(jj = 0; jj < this->n_dofs; ++jj){
			this->k_e(ii,jj) = 0.0;
		}
	}

}

void IntegrationHelper::outElement()
{
	/*perform clean up and destroy specific field element*/
	apf::destroyElement(this->field_element);
}

void IntegrationHelper::atPoint(apf::Vector3 const& p, double w, double dV)
{
	/**/
	std::cout << "p: " << p << std::endl;
	std::cout << "weight: " << w << " dV: " << dV << std::endl;

	apf::MeshElement* me = apf::getMeshElement(this->field_element);
	apf::Vector3 x;
	apf::mapLocalToGlobal(me, p, x);
	std::cout << "global x: " << x << std::endl;
}
