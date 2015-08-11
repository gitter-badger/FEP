#include <iostream>

#include <apfMesh.h>

#include "TestClass.h"

TestClass::TestClass(int order)
{
	this->secret = order;
}


ChildClass::ChildClass(int order, int cats) : TestClass(order)
{
	this->key  = cats;
}

void ChildClass::someMethod()
{
	std::cout << "child method" << std::endl;
}

Integrate::Integrate(uint32_t order, apf::Field* field) :
  apf::Integrator(order),
  u(field)
{
  ndims = apf::getMesh(u)->getDimension();
}

void Integrate::inElement(apf::MeshElement* me)
{
  e = apf::createElement(u,me);
  ndofs = apf::countNodes(e);
  fe.setSize(ndofs);
  ke.setSize(ndofs,ndofs);
  for (int a=0; a < ndofs; ++a)
  {
    fe(a) = 0.0;
    for (int b=0; b < ndofs; ++b)
      ke(a,b) = 0.0;
  }
}

void Integrate::outElement()
{
  apf::destroyElement(e);
}

void Integrate::atPoint(apf::Vector3 const& p, double w, double dV)
{
  std::cout << "p: " << p << std::endl;
  std::cout << "weight: " << w << " dV: " << dV << std::endl;

  apf::MeshElement* me = apf::getMeshElement(this->e);
  apf::Vector3 x;
  apf::mapLocalToGlobal(me, p, x);
  std::cout << "global x: " << x << std::endl;
  // apf::NewArray<double> BF;
  // apf::getShapeValues(e,p,BF);

  // apf::NewArray<apf::Vector3> gradBF;
  // apf::getShapeGrads(e,p,gradBF);

  // apf::Vector3 x;
  // apf::MeshElement* me = apf::getMeshElement(e);
  // apf::mapLocalToGlobal(me,p,x);

  // for (int a=0; a < ndofs; ++a)
  // {
  //   fe(a) += rhs(x) * BF[a] * w * dv;
  //   for (int b=0; b < ndofs; ++b)
  //   for (int i=0; i < ndims; ++i)
  //     ke(a,b) += gradBF[a][i] * gradBF[b][i] * w * dv; 
  // }
}
