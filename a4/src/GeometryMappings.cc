#include "GeometryMappings.h"

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
	void(*fnc_ptr)(apf::MeshEntity* e,
					std::vector< uint64_t > & nodes,
					std::vector < double > & d))
{

}

