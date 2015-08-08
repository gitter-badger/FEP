#include "ElasticAnalysis2D.h"

ElasticAnalysis2D::ElasticAnalysis2D() {

}

ElasticAnalysis2D::~ElasticAnalysis2D() {

}

uint32_t ElasticAnalysis2D::setup() {
	return 0;
}

uint32_t ElasticAnalysis2D::solve() {
	return 0;
}

uint32_t ElasticAnalysis2D::makeStiffnessContributor(apf::MeshEntity* e) {
	return 0;
}

uint32_t ElasticAnalysis2D::makeForceContributor(apf::MeshEntity* e) {
	return 0;
}

uint32_t ElasticAnalysis2D::makeConstraint(apf::MeshEntity* e) {
	return 0;
}

uint32_t ElasticAnalysis2D::recover() {
	return 0;
}
