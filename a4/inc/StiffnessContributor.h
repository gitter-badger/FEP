#ifndef STIFFNESS_CONTRIBUTOR_H
#define STIFFNESS_CONTRIBUTOR_H 1

#include <stdint.h>
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>

class StiffnessContributor : public apf::Integrator
{
public:
	StiffnessContributor(apf::Field *f, uint32_t integrate_order);
	void inElement(apf::MeshElement *me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dV);
    apf::DynamicVector fe;
    apf::DynamicMatrix ke;
private:
    uint32_t ndofs;
    uint32_t ndims;
    apf::Field* field;
    apf::Element* field_element;

};

#endif
