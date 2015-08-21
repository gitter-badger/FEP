#ifndef TEST_MESH_BUILDER
#define TEST_MESH_BUILDER
#include <stdint.h>

#include <apf.h>
#include <apfMesh2.h>

/*this numbering is left on the mesh for the purpose of unittesting
* the edge boundary conditions end up where they are supposed to
* this is for DEVELOPEMENT ONLY*/
#define SECRET_BUILDER_NUMBERING "SecretBuilderNumbering"

class MeshBuilder
{
public:
	MeshBuilder();
	~MeshBuilder();

	void build2DRectQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
		uint32_t y_elms, double x0, double y0, double xf, double yf);
	void build2DJaggedQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
		uint32_t y_elms, double x0, double y0, double xf, double yf);

	void build2DRectTriMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
		uint32_t y_elms, double x0, double y0, double xf, double yf);
	
	void AdjReorder(apf::Mesh2* & mesh);

private:

};

class EntityNumberer : public apf::BuildCallback
{
public:
	EntityNumberer(apf::Mesh* mesh);
	~EntityNumberer();

	void call(apf::MeshEntity *e);

	/*one must configure these two below before using the call back object
	* this is a hacky way to make sure we numbe the one we want*/
	int number_to_apply;	
	apf::MeshTag* tag;
private:
	apf::Mesh* mesh;

};

#endif
