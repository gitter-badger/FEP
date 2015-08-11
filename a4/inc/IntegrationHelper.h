#ifndef INTEGRATION_HELPER
#define INTEGRATION_HELPER

#include <stdint.h>

#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>

class IntegrationHelper : public apf::Integrator
{
	/*This class caches a user defined function to be integrated over
	* the domain of a mesh element. The fuction is passed in through 
	* the constructor. This compute the element stiffness matrix
	*/
  public:
    IntegrationHelper(apf::Field* f, uint32_t order);
    void inElement(apf::MeshElement*);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dv);
    apf::DynamicVector f_e;
    apf::DynamicMatrix k_e;
  private:
    int n_dofs;
    int n_dimensions;
    apf::Field* field;
    apf::Element* field_element;
};

#endif
