#ifndef MESH_ADJ_REODRDER
#define MESH_ADJ_REODRDER 

#include <apfNumbering.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <stdint.h>

void adjReorder(apf::Mesh* m, apf::FieldShape* shape, uint32_t components,
            apf::Numbering* & inNodes, apf::Numbering* & elmNodes);
#endif
