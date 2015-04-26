#ifndef TEST_MESH_BUILDER
#define TEST_MESH_BUILDER
#include <stdint.h>

#include <apfMesh2.h>


class MeshBuilder
{
public:
	MeshBuilder();
	~MeshBuilder();

	void build2DRectQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
		uint32_t y_elms, double x0, double y0, double xf, double yf);
	void build2DJaggedQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
		uint32_t y_elms, double x0, double y0, double xf, double yf);

	void build2DTriQuadMesh(const apf::Mesh2* mesh, uint32_t xunits, 
		uint32_t yunits, double x0, double y0, double xf, double yf);
	void AdjReorder(apf::Mesh2* & mesh);

private:

};

#endif
