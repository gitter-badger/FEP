#include "GeometryMappings.h"

void noConstraint(
	apf::MeshElement *e,
	apf::Numbering* nodeNums,
	std::vector<uint64_t> &dofs,
	std::vector<double> & disp)
{
	/*empty the constraint vector*/
	dofs.clear();
	disp.clear();
}

void zeroDisplacementX_2D(
	apf::MeshElement *me,
	apf::Numbering* nodeNums,
	std::vector<uint64_t> & dofs,
	std::vector<double> & disp)
{

	int entity_type = me->getMesh()->getType(me->getEntity());
	uint32_t nnodes = apf::countElementNodes(me->getMesh()->getShape(), entity_type);
	/*in 2D we assume there are only 2 dofs per node*/
	apf::NewArray< int > node_mapping(nnodes*2);
	uint32_t tmp_sz = apf::getElementNumbers(nodeNums, me->getEntity(), node_mapping);
	assert((nnodes*2) == tmp_sz);
	/*resize the vector to hold only the fixed dofs in this case is exactly
	* nnodes*/
	dofs.resize(nnodes);
	disp.resize(nnodes);
	std::size_t curs_indx = 0;
	for(std::size_t ii = 0; ii < tmp_sz; ++ii) {
		/*we assume the ordering of local dofs is X_0, Y_0, X_1, Y_1 ...
		* therefor X terms are the even terms*/
		if(0 == ii %2) {
			dofs[curs_indx] = node_mapping[ii];
			disp[curs_indx] = 0.0;
			curs_indx++;
		}
	}
	assert(curs_indx == nnodes);
}

void zeroDisplacementY_2D(
	apf::MeshElement *me,
	apf::Numbering* nodeNums,
	std::vector<uint64_t> & dofs,
	std::vector<double> & disp)
{
	int entity_type = me->getMesh()->getType(me->getEntity());
	uint32_t nnodes = apf::countElementNodes(me->getMesh()->getShape(), entity_type);
	/*in 2D we assume there are only 2 dofs per node*/
	apf::NewArray< int > node_mapping(nnodes*2);
	uint32_t tmp_sz = apf::getElementNumbers(nodeNums, me->getEntity(), node_mapping);
	assert((nnodes*2) == tmp_sz);
	/*resize the vector to hold only the fixed dofs in this case is exactly
	* nnodes*/
	dofs.resize(nnodes);
	disp.resize(nnodes);
	std::size_t curs_indx = 0;
	for(std::size_t ii = 0; ii < tmp_sz; ++ii) {
		/*we assume the ordering of local dofs is X_0, Y_0, X_1, Y_1 ...
		* therefor y terms are the odd terms*/
		if(1 == ii %2) {
			dofs[curs_indx] = node_mapping[ii];
			disp[curs_indx] = 0.0;
			curs_indx++;
		}
	}
	assert(curs_indx == nnodes);
}

GeometryMappings::GeometryMappings()
{
	nuemann_keys.clear();
	dirchelet_keys.clear();
	neumann_map.clear();
	dirchelet_map.clear();
	/*we insert some simple boundary conditions*/
}

GeometryMappings::~GeometryMappings()
{

}

void GeometryMappings::addNeumannMapping(
	uint64_t key,
	apf::Vector3(*neumann_fnc)(apf::Vector3 const& p))
{

}

void GeometryMappings::addDircheletMapping(
	uint64_t key,
	void(*fnc_ptr)(apf::MeshElement* me,
					apf::Numbering* nodeNums,
					std::vector< uint64_t > & nodes,
					std::vector < double > & d))
{

}

