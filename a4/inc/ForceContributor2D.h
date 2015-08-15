#ifndef STIFFNESS_CONTRIBUTOR_H
#define STIFFNESS_CONTRIBUTOR_H 1

#include <vector>
#include <stdint.h>
#include <apf.h>


class ForceContributor2D : public apf::Integrator
{
public:
	ForceContributor2D(apf::Field *f, uint32_t integrate_order, apf::Vector3(*fnc)(apf::Vector3 const& p));
	void inElement(apf::MeshElement *me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dV);
    std::vector<double> fe;
private:
    uint32_t ndofs;
    uint32_t nnodes;
    uint32_t ndims;
    apf::Field* field;
    apf::Element* field_element;
    apf::Vector3(*fnc)(apf::Vector3 const& p);
};

#endif
