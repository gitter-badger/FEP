#ifndef MESH_ADJ_REODRDER
#define MESH_ADJ_REODRDER 

#include <apfNumbering.h>
#include <apfMesh.h>
#include <apfShape.h>

void adjReorder(apf::Mesh* m, apf::FieldShape* shape, int components,
            apf::Numbering* & inNodes, apf::Numbering* & elmNodes);
#endif
