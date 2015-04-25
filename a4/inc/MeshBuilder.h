#ifndef TEST_MESH_BUILDER
#define TEST_MESH_BUILDER
#include <stdint.h>

#include <apfMesh2.h>


class MeshBuilder
{
public:
	MeshBuilder();
	~MeshBuilder();

	void build2DRectQuadMesh(const apf::Mesh2* mesh, uint32_t xunits, 
		uint32_t yunits, uint32_t x0, uint32_t y0, uint32_t xf, uint32_t vf);
	void build2DTriQuadMesh(const apf::Mesh2* mesh, uint32_t xunits, 
		uint32_t yunits, uint32_t x0, uint32_t y0, uint32_t xf, uint32_t vf);

private:

};

#endif
