#ifndef TEST_CLASS_H
#define TEST_CLASS_H

#include <apf.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <stdint.h>

class TestClass
{
public:
	TestClass(int order);
	~TestClass();
	virtual void someMethod() = 0;

	int secret;
};

class ChildClass : public TestClass
{
public:
	ChildClass(int order, int cats);
	~ChildClass();

	void someMethod();
	/* data */
	int key;
};


struct IntegrateInput
{
  uint32_t order;
  apf::Field* field;
};

class Integrate : public apf::Integrator
{
  public:
    Integrate(uint32_t order, apf::Field* field);
    void inElement(apf::MeshElement*);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dv);
    apf::DynamicVector fe;
    apf::DynamicMatrix ke;
  private:
    int ndofs;
    int ndims;
    apf::Field* u;
    apf::Element* e;
};


#endif
