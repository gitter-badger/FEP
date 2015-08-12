#ifndef STIFFNESS_CONTRIBUTOR_H
#define STIFFNESS_CONTRIBUTOR_H 1

#include <stdint.h>
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfDynamicVector.h>

class StiffnessContributor2D : public apf::Integrator
{
public:
	StiffnessContributor2D(apf::Field *f, apf::Matrix< 3, 3 > D, uint32_t integrate_order);
	void inElement(apf::MeshElement *me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dV);
    void computeBtranspose(apf::Vector3 const& grad, apf::Matrix< 2,3 > & B_T);
    void computeB(apf::Vector3 const& grad, apf::Matrix< 3,2 > & B);
    apf::DynamicVector fe;
    apf::DynamicMatrix ke;
private:
    uint32_t ndofs;
    uint32_t nnodes;
    uint32_t ndims;
    apf::Field* field;
    apf::Element* field_element;
    apf::Matrix< 3,3 > D;
};

#endif
