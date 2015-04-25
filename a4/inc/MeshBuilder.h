#ifndef TEST_MESH_BUILDER
#define TEST_MESH_BUILDER
#include <apfMesh2.h>

class TestMeshBuilder
{
public:
	TestMeshBuilder();
	~TestMeshBuilder();

	void build2DRectQuadMesh(const apf::Mesh2* mesh, uint32_t xunits, 
		uint32_t yunits, uint32_t x0, uint32_t y0, uint32_t xf, uint32_t vf);
	void build2DTriQuadMesh(const apf::Mesh2* mesh, uint32_t xunits, 
		uint32_t yunits, uint32_t x0, uint32_t y0, uint32_t xf, uint32_t vf);

private:

};

#endif
