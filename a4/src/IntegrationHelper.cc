#include "IntegrationHelper.h"

IntegrationHelper::IntegrationHelper(apf::Field* field, double(*func)(apf::Vector3 const& p), uint32_t order) : apf::Integrator(order)
{
	this->field = field;
	this->fnc_ptr = func;
	this->n_dimensions = apf::getMesh(field)->getDimension();
}

IntegrationHelper::IntegrationHelper(apf::Field* field, apf::Vector3(*func)(apf::Vector3 const& p), uint32_t order) : apf::Integrator(order)
{
	this->field = field;
	this->fnc_ptr = func;
	this->n_dimensions = apf::getMesh(field)->getDimension();
}

void IntegrationHelper::inElement(apf::MeshElement* me)
{
	/*create the Field Element from the MeshElement*/
	this->field_element = apf::createElement(this->field, me);
	/*determine the size of the force matrices*/
	this->n_dofs = apf::countNodes(this->field_element);
	this->f_e.setSize(this->n_dofs);
	/*zero the internal matricies*/
	uint32_t ii, jj;
	for(ii = 0; ii < this->n_dofs; ++ii) {
		this->f_e = 0.0;
		for(jj = 0; jj < this->n_dofs; ++ii){
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
	/*Purpose: Compute force contributions for function*/
	/*storage for*/
	apf::NewArray<double> 

	/*location of point p in global coordinates*/
	apf::Vector3 global_p;

}
