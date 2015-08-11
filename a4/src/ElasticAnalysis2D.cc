#include "ElasticAnalysis2D.h"

ElasticAnalysis2D::ElasticAnalysis2D(apf::Mesh* m)
{
	this->m = m;
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

	if(entity_type == apf::Mesh::QUAD || entity_type == apf::Mesh::TRIANGLE){
		std::cout << "face" << std::endl;
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
