#include <stdint.h>

#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
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

void MeshBuilder::build2DRectQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
	uint32_t y_elms, double x0, double y0, double xf, double yf)
{
	/*there is nothing to be done for a zero size mesh so return*/
	/*there is nothing to be done for a zero size mesh so return*/
	if ((x_elms == 0) || (y_elms == 0)) {
		mesh = NULL;
		return;
	}
	gmi_register_null();
  	gmi_model* g = gmi_load(".null");
  	mesh = apf::makeEmptyMdsMesh(g, 2, false);

  	uint64_t total_elms = (x_elms + 1) * (y_elms + 1);
  	long double x_size = xf - x0;
  	long double y_size = yf - y0;

  	//create an empty numbering
	apf::Numbering* numbers = apf::createNumbering(mesh,"my_numbers", 
												   mesh->getShape(), 1);
	int node_number = 0;
	//create an array to hold all our vertices
	apf::MeshEntity** vertices = new apf::MeshEntity*[total_elms];
	for(uint32_t counter = 0; counter < total_elms; counter++) {
		vertices[counter] = mesh->createVert(0);
	}
	//create a pointer to pass in quad vertices
	apf::MeshEntity* quad_verts[4];

	//rather than test for that zero,zero element every time
	//just instantiate it manually once
	apf::Vector3* temp_vec = new apf::Vector3();
	double xyz[3] = {x0, y0, 0};
	temp_vec->fromArray(xyz);
	int vert_index = 0; // use this to place verts in right spot
	mesh->setPoint(vertices[vert_index],0,*(temp_vec));
	//numbering each node in order of creation
	apf::number(numbers,vertices[vert_index],0,0, ++node_number);
	
	for(uint32_t y_c = 1; y_c <= y_elms; ++y_c) {
		//for the first column of x we make element to left
		xyz[0] = x0;
		xyz[1] = (static_cast<double>(y_c) * y_size / 
				  static_cast<double> (y_elms)) + y0;
		temp_vec->fromArray(xyz); //load the values
		vert_index = (x_elms + 1) * y_c;
		//print the location of each vertex
		mesh->setPoint(vertices[vert_index],0,*(temp_vec));
		apf::number(numbers,vertices[vert_index],0,0, ++node_number);

		for(uint32_t x_c = 1; x_c <= x_elms; x_c++) {
			xyz[0] = (static_cast<double>(x_c) * x_size / 
					  static_cast<double>(x_elms)) + x0;
			//for the first row of y we construct an extra element below
			if(y_c == 1) {
				xyz[1] = y0;
				temp_vec->fromArray(xyz);
				vert_index = x_c;
				mesh->setPoint(vertices[vert_index],0,*(temp_vec));
				apf::number(numbers,vertices[vert_index],0,0,++node_number);
			} //now create the actual row
			xyz[1] = (static_cast<double>(y_c) * y_size / 
					  static_cast<double>(y_elms)) + y0;
			temp_vec->fromArray(xyz);
			vert_index = y_c*(x_elms+1) + x_c;
			mesh->setPoint(vertices[vert_index],0,*(temp_vec));
			apf::number(numbers,vertices[vert_index],0,0,++node_number);
			//print the location of the vertex
			//create the quad element
			quad_verts[0] = vertices[vert_index];
			quad_verts[1] = vertices[(vert_index - 1)];
			quad_verts[2] = vertices[(vert_index - x_elms - 2)];
			quad_verts[3] = vertices[(vert_index - x_elms - 1)];
			apf::buildElement(mesh, 0, apf::Mesh::QUAD, quad_verts);
		}
	}
	apf::deriveMdsModel(mesh);//this makes CAD model for classification
	//accept the changes
	mesh->acceptChanges();
	mesh->verify();
	delete[] vertices;

}
void build2DJaggedQuadMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
		uint32_t y_elms, double x0, double y0, double xf, double yf) 
{


}
void MeshBuilder::build2DRectTriMesh(apf::Mesh2* & mesh, uint32_t x_elms, 
	uint32_t y_elms, double x0, double y0, double xf, double yf) 
{
	/*there is nothing to be done for a zero size mesh so return*/
	if ((x_elms == 0) || (y_elms == 0)) {
		mesh = NULL;
		return;
	}
	gmi_register_null();
  	gmi_model* g = gmi_load(".null");
  	mesh = apf::makeEmptyMdsMesh(g, 2, false);

  	uint64_t total_elms = (x_elms + 1) * (y_elms + 1);
  	long double x_size = xf - x0;
  	long double y_size = yf - y0;

  	//create an empty numbering
	apf::Numbering* numbers = apf::createNumbering(mesh,"my_numbers", 
												   mesh->getShape(), 1);
	int node_number = 0;
	//create an array to hold all our vertices
	apf::MeshEntity** vertices = new apf::MeshEntity*[total_elms];
	for(uint32_t counter = 0; counter < total_elms; counter++) {
		vertices[counter] = mesh->createVert(0);
	}
	//create a pointer to pass in quad vertices
	apf::MeshEntity* tri_verts[3];

	//rather than test for that zero,zero element every time
	//just instantiate it manually once
	apf::Vector3* temp_vec = new apf::Vector3();
	double xyz[3] = {x0, y0, 0};
	temp_vec->fromArray(xyz);
	int vert_index = 0; // use this to place verts in right spot
	mesh->setPoint(vertices[vert_index],0,*(temp_vec));
	//numbering each node in order of creation
	apf::number(numbers,vertices[vert_index],0,0, ++node_number);
	
	for(uint32_t y_c = 1; y_c <= y_elms; ++y_c) {
		//for the first column of x we make element to left
		xyz[0] = x0;
		xyz[1] = (static_cast<double>(y_c) * y_size / 
				  static_cast<double> (y_elms)) + y0;
		temp_vec->fromArray(xyz); //load the values
		vert_index = (x_elms + 1) * y_c;
		//print the location of each vertex
		mesh->setPoint(vertices[vert_index],0,*(temp_vec));
		apf::number(numbers,vertices[vert_index],0,0, ++node_number);

		for(uint32_t x_c = 1; x_c <= x_elms; x_c++) {
			xyz[0] = (static_cast<double>(x_c) * x_size / 
					  static_cast<double>(x_elms)) + x0;
			//for the first row of y we construct an extra element below
			if(y_c == 1) {
				xyz[1] = y0;
				temp_vec->fromArray(xyz);
				vert_index = x_c;
				mesh->setPoint(vertices[vert_index],0,*(temp_vec));
				apf::number(numbers,vertices[vert_index],0,0,++node_number);
			} //now create the actual row
			xyz[1] = (static_cast<double>(y_c) * y_size / 
					  static_cast<double>(y_elms)) + y0;
			temp_vec->fromArray(xyz);
			vert_index = y_c*(x_elms+1) + x_c;
			mesh->setPoint(vertices[vert_index],0,*(temp_vec));
			apf::number(numbers,vertices[vert_index],0,0,++node_number);
			//print the location of the vertex
			//creat two triangular elements from quad
			tri_verts[0] = vertices[vert_index];
			tri_verts[1] = vertices[(vert_index - 1)];
			tri_verts[2] = vertices[(vert_index - x_elms - 2)];
			apf::buildElement(mesh, 0, apf::Mesh::TRIANGLE, tri_verts);
			tri_verts[0] = vertices[(vert_index - x_elms - 2)];
			tri_verts[1] = vertices[(vert_index - x_elms - 1)];
			tri_verts[2] = vertices[vert_index];
			apf::buildElement(mesh, 0, apf::Mesh::TRIANGLE, tri_verts);
		}
	}
	apf::deriveMdsModel(mesh);//this makes CAD model for classification
	//accept the changes
	mesh->acceptChanges();
	mesh->verify();
	delete[] vertices;

}

void AdjReorder(apf::Mesh2* & mesh) {

}
