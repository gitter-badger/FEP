#include "ElasticAnalysis2D.h"
#include "StiffnessContributor2D.h"
#include "ForceContributor2D.h"
#include "MeshAdjReorder.h"

#include <fstream>

#include <PCU.h>

#define NUM_COMPONENTS 2

apf::Vector3 dummy(apf::Vector3 const& p){
	return apf::Vector3(10,0,0);
}

apf::Matrix< 3,3 > buildD(double E, double Nu, bool use_plane_stress) {
	apf::Matrix< 3,3 > D;
	/*find the Lame parameters for plane strain*/
	double lambda, mu;
	lambda = (Nu * E) / ((1.0 + Nu) * (1.0 - 2.0 * Nu));
	mu = E / (2.0 + 2.0 * Nu);
	/*modify for plane stress*/
	if(use_plane_stress) {
		lambda = (2.0 * lambda * mu) / ( lambda + 2.0 * mu);
	}
	/*initialize the D matrix for which ever case we are using*/
	D[0][0] = lambda + 2 * mu;
	D[0][1] = lambda;
	D[0][2] = 0.0;
	D[1][0] = lambda;
	D[1][1] = lambda + 2 * mu;
	D[1][2] = 0.0;
	D[2][0] = 0.0;
	D[2][1] = 0.0;
	D[2][2] = mu;
	return D;
}

ElasticAnalysis2D::ElasticAnalysis2D(struct ElasticAnalysisInput & in) :
	integration_order(in.integration_order),
	m(in.m)
{
	/*we create a 3 dimensional field but only will use 2 dimensions of that, (x & y)*/
	this->field = createField(this->m, SOLUTION_FIELD_NAME, apf::VECTOR, this->m->getShape());
	apf::zeroField(this->field);

	bool use_plane_stress = true;
	this->D = buildD(in.E, in.Nu, use_plane_stress);

	/*set up element numberings we will use for assembly*/
	this->nodeNums = apf::createNumbering(this->m, NODE_NUM_TAG_NAME, this->m->getShape(), NUM_COMPONENTS);
    this->faceNums = apf::createNumbering(this->m, FACE_NUM_TAG_NAME, apf::getConstant(this->m->getDimension()), 1);
	if(in.reorder == true){
		adjReorder(this->m, this->m->getShape(), NUM_COMPONENTS, this->nodeNums, this->faceNums);
	}
	/*compute the global degrees of freedom for the mesh*/
	std::size_t n_global_dofs = apf::countNodes(this->nodeNums) * NUM_COMPONENTS;
	/*initialize the linear system*/
	this->linsys = new AlgebraicSystem(n_global_dofs);
}

ElasticAnalysis2D::~ElasticAnalysis2D()
{
	apf::destroyNumbering(this->nodeNums);
	apf::destroyNumbering(this->faceNums);
	delete this->linsys;
}

uint32_t ElasticAnalysis2D::setup()
{
	apf::MeshIterator* it;
	apf::MeshEntity* e;
	/*find all of the boundary conditions fixed,
	* chech all edges classified on model edge and
	* all vertices classified on model vertex*/
	bool _arb_flag = true;

	it = this->m->begin(1);
	while((e = this->m->iterate(it))) {
		if(_arb_flag){ 
			// this->makeConstraint(e);
			_arb_flag = false;
		}
	}
	/*randomly fix only the last edge*/

	this->m->end(it);
	it = this->m->begin(0);
	while((e = this->m->iterate(it))) {

	}
	this->m->end(it);

	/*now that all fixed dofs have been accounted for,
	* allow assembly of force and stiffness contributors*/
	this->linsys->beginAssembly();

	/*iterate over the faces*/
	it = this->m->begin(2);
	while((e = this->m->iterate(it))) {
		this->makeStiffnessContributor(e);
		this->makeForceContributor(e);
	}
	this->m->end(it);
	/*then pick up the edges*/
	it = this->m->begin(1);
	while((e = this->m->iterate(it))) {
		this->makeForceContributor(e);
	}
	this->m->end(it);

	/*now synchronize the linear system in preparation for solving*/
	this->linsys->synchronize();
	return 0;
}

uint32_t ElasticAnalysis2D::solve()
{
	std::fstream myfile;
	myfile.open("K.txt");
	myfile << this->linsys->K;
	myfile.close();

	double t0 = PCU_Time();
	this->linsys->solve();
	double t1 = PCU_Time();
	std::cout << "System solved in: " << t1- t0 << " seconds" <<std::endl;
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
		apf::destroyMeshElement(me);
		/*view the intermediate matrix*/
		uint32_t n_l_dofs = apf::countElementNodes(this->m->getShape(), entity_type) * NUM_COMPONENTS;
		apf::NewArray< int > node_mapping(n_l_dofs);
		apf::getElementNumbers(nodeNums, e, node_mapping);
		// std::cout  << "length of numbering: " << length << std::endl;
		// for(uint32_t ii = 0; ii < nnodes; ++ii){
		// 	std::cout << "Node " << ii << ": " << node_mapping[ii] << std::endl;
		// }
		this->linsys->assemble(stiff.ke, node_mapping, n_l_dofs);

	} else {
		/*only accepts faces, so indicate improper input*/
		std::cout << "entity type: " << entity_type << std::endl;
		return 1;
	}
	return 0;
}

uint32_t ElasticAnalysis2D::makeForceContributor(apf::MeshEntity* e)
{
	int entity_type = this->m->getType(e);
	apf::Vector3(*fnc_ptr)(apf::Vector3 const& p);
	fnc_ptr = &dummy;
	ForceContributor2D force(this->field, this->integration_order, fnc_ptr);

	apf::MeshElement* me = apf::createMeshElement(this->m, e);
	force.process(me);
	apf::destroyMeshElement(me);

	/*view the intermediate matrix*/
	uint32_t n_l_dofs = apf::countElementNodes(this->m->getShape(), entity_type) * NUM_COMPONENTS;
	apf::NewArray< int > node_mapping(n_l_dofs);
	apf::getElementNumbers(nodeNums, e, node_mapping);

	this->linsys->assemble(force.fe, node_mapping, n_l_dofs);

	return 0;
}

uint32_t ElasticAnalysis2D::makeConstraint(apf::MeshEntity* e)
{
	int entity_type = this->m->getType(e);
	uint32_t n_l_dofs = apf::countElementNodes(this->m->getShape(), entity_type) * NUM_COMPONENTS;
	apf::NewArray< int > node_mapping(n_l_dofs);
	apf::getElementNumbers(nodeNums, e, node_mapping);

	std::vector< uint32_t > fixed_mapping(0);
	std::vector< double > displacement(0);


	/*right now all fixed nodes will have zero displacement in the x*/
	for(uint32_t ii = 0; ii < n_l_dofs; ++ii) {
		if(node_mapping[ii] % 2 == 0) {
			fixed_mapping.push_back(node_mapping[ii]);
			displacement.push_back(0.0);
		}
	}
	this->linsys->addBoundaryConstraint(displacement, fixed_mapping);

	return 0;
}

uint32_t ElasticAnalysis2D::recover()
{
	this->displacement.clear();
	this->strain.clear();
	this->stress.clear();

	this->linsys->extractDisplacement(this->displacement);
	/*resize each elements and default construct that */

	return 0;
}
