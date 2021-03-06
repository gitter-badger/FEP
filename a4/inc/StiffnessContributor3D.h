#ifndef STIFFNESS_CONTRIBUTOR_H
#define STIFFNESS_CONTRIBUTOR_H 1

#include <stdint.h>
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>
#include <apfMatrix.h>

class StiffnessContributor3D : public apf::Integrator
{
public:
	StiffnessContributor3D(apf::Field *f, apf::Matrix< 6,6 > D, uint32_t integrate_order);
	virtual void inElement(apf::MeshElement *me);
    virtual void outElement();
    virtual void atPoint(apf::Vector3 const& p, double w, double dV);
    apf::DynamicVector fe;
    apf::DynamicMatrix ke;
private:
    uint32_t ndofs;
    uint32_t ndims;
    apf::Field* field;
    apf::Element* field_element;
    apf::Matrix< 6,6 > D;
};

#endif
