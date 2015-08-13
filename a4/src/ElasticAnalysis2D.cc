#include "ElasticAnalysis2D.h"
#include "StiffnessContributor2D.h"

ElasticAnalysis2D::ElasticAnalysis2D(
	apf::Mesh* m,
	uint32_t integration_order,
	double E,
	double Nu)
{
	this->integration_order = integration_order;
	this->m = m;
	this->field = createField(this->m, "Field_1", apf::VECTOR, this->m->getShape());
	apf::zeroField(this->field);
	/*find the Lame parameters for plane strain*/
	double lambda, mu;
	lambda = (Nu * E) / ((1.0 + Nu) * (1.0 - 2.0 * Nu));
	mu = E / (2.0 + 2.0 * Nu);
	/*modify for plane stress*/
	const uint32_t PLANE_STRESS = 1;
	if(PLANE_STRESS) {
		lambda = (2.0 * lambda * mu) / ( lambda + 2.0 * mu);
	}
	/*initialize the D matrix for which ever case we are using*/
	this->D[0][0] = lambda + 2 * mu;
	this->D[0][1] = lambda;
	this->D[0][2] = 0.0;
	this->D[1][0] = lambda;
	this->D[1][1] = lambda + 2 * mu;
	this->D[1][2] = 0.0;
	this->D[2][0] = 0.0;
	this->D[2][1] = 0.0;
	this->D[2][2] = mu;
}

ElasticAnalysis2D::~ElasticAnalysis2D()
{

}

uint32_t ElasticAnalysis2D::setup()
{
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	/*iterate over the faces first*/
	it = this->m->begin(2);
	while((e = this->m->iterate(it))){
		this->makeStiffnessContributor(e);
	}
	return 0;
}

uint32_t ElasticAnalysis2D::solve()
{
	return 0;
}

uint32_t ElasticAnalysis2D::makeStiffnessContributor(apf::MeshEntity* e)
{

	int entity_type = this->m->getType(e);
	/*stiffness contributors are entirely local to this method, and are
	* assembled directly as they are generated. For meshes of same element
	* type the allocation will not change from element to element. For mixed
	* meshes, the interface will remain the same, but the allocation size will
	* be stuck on the largest element size. For example if the largest element
	* is a 4th order Lagrange element then the size of the local stiffness
	* contributor will grow as the entities are processed up to the size of
	* the largest element, and after that it will remain this size until the
	* end of the program. This is acceptable right now since we do
	* intend to do mixed meshes with different orders. Even if they are
	* the same order, the difference in size between different elements
	* will increase as the order of the element shape functions increase.
	* We will stick to low order (1 or 2) shape functions so this size
	* descrapancy is minimal. For comparison see table below for number of
	* lagrange nodes of Nth order for triangular and quadrilateral elements
	* order 	1 	2 	3	4
	* tri  		3 	6   10  15
	* quad 		4   9	16	25
	* As one can see. the difference between a triangular and quadrilateral
	* element is a difference of ((50 dofs)^2 - (30 dofs)^2) is allocating
	* 78% more storage than is actually required for processing the triangular
	* element. Clearly a strong canidate for causing unneccesary cache
	* misses.
	**/
	StiffnessContributor2D stiff(this->field, this->D, this->integration_order);

	if(entity_type == apf::Mesh::QUAD || entity_type == apf::Mesh::TRIANGLE){
		apf::MeshElement* me = apf::createMeshElement(this->m, e);
		stiff.process(me);

		/*view the intermediate matrix*/
		std::cout << stiff.ke << std::endl << "=====================" << std::endl;

	} else {
		/*only accepts faces, so indicate improper input*/
		std::cout << "entity type: " << entity_type << std::endl;
		return 1;
	}
	return 0;
}

uint32_t ElasticAnalysis2D::makeForceContributor(apf::MeshEntity* e)
{
	return 0;
}

uint32_t ElasticAnalysis2D::makeConstraint(apf::MeshEntity* e)
{
	return 0;
}

uint32_t ElasticAnalysis2D::recover()
{
	return 0;
}
