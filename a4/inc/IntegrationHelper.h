#ifndef INTEGRATION_HELPER
#define INTEGRATION_HELPER 1

#include <stdint.h>

#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>

class IntegrationHelper : public apf::Integrator
{
	/*This class caches a user defined function to be integrated over
	* the domain of a mesh element. The fuction is passed in through 
	* the constructor. This compute the element stiffness matrix
	**/
public:
	/*scalar valued functions*/
	IntegrationHelper(apf::Field* field, double(*func)(apf::Vector3 const& p), uint32_t order);
	/*vector valued functions*/
	IntegrationHelper(apf::Field* field, apf::Vector3(*func)(apf::Vector3 const& p), uint32_t order);
	~IntegrationHelper();

	void inElement(apf::MeshElement* me);
	void outElement();
	void atPoint(apf::Vector3 const&p, double w, double dV);
	apf::apfDynamicVector f_e;

private:
	apf::Field *field;
	apf::Element* field_element;
	uint32_t n_dimensions;
	uint32_t n_dofs;
	/*fuction pointer we will evalutate, overloaded to support different
	* type of functions*/
	double(*fnc_ptr)(apf::Vector3 const& p);
	apf::Vector3(*fnc_ptr)(apf::Vector3 const& p);
};

#endif
