#include <stdint.h>

#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include "MeshBuilder.h"

MeshBuilder::MeshBuilder() {
}
MeshBuilder::~MeshBuilder() {
}

void MeshBuilder::build2DRectQuadMesh(const apf::Mesh2* mesh, uint32_t xunits, 
	uint32_t yunits, uint32_t x0, uint32_t y0, uint32_t xf, uint32_t vf) {
	printf("Inside quad mesh builder\n");
}
void MeshBuilder::build2DTriQuadMesh(const apf::Mesh2* mesh, uint32_t xunits, 
	uint32_t yunits, uint32_t x0, uint32_t y0, uint32_t xf, uint32_t vf) {
	printf("Inside triagle mesh builder\n");
}
