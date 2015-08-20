#ifndef GEOMETRY_MAPPING_H
#define GEOMETRY_MAPPING_H 1

#include <map>
#include <vector>
#include <set>
#include <stdint.h>

#include <apf.h>
#include <apfMesh.h>

/*include some simple defaults*/
enum GeoMapOptions {
	FIXED_X = 0,
	FIXED_Y = 1,
	FIXED_XY = 1
};

typedef apf::Vector3(*neumann_fnc)(apf::Vector3 const& );
typedef void(*bound_gen)(apf::MeshEntity* , std::vector< uint64_t >&, std::vector < double > & );

class GeometryMappings
{
	/*manually configured in the .cc file, this applies
	* all boundary conditions to geometric model based on integer tags
	*/

public:
	GeometryMappings();
	~GeometryMappings();

	void addNeumannMapping(uint64_t key, apf::Vector3(*neumann_fnc)(apf::Vector3 const& ));
	void addDircheletMapping(uint64_t key, void(*fnc_ptr)(apf::MeshEntity* , std::vector< uint64_t > & nodes, std::vector < double > & d));

	std::map< uint64_t,neumann_fnc > neumann_map;
	std::map< uint64_t,bound_gen > dirchelet_map;

private:
	std::set< uint64_t > nuemann_keys;
	std::set< uint64_t > dirchelet_keys;

};


#endif

