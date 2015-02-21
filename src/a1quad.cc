#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>

#include <stdio.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 2, false);
  /*Create a uniform 2D quadirlateral mesh of 4 elements by 6 elements on a
    square of edge length 8.0, with lower left vertex placed at 0,0
    provide a list that gives the coordinates for each vertex. Provide
    a list of the vertex numbers for each element. Provide a copy of an
    image of the mesh
   */
  int x_elms = 4;
  int y_elms = 6;
  int total_elms = (x_elms+1) * (y_elms+1);
  double x_size = 8.0;
  double y_size = 8.0;
  //create an empty numbering
  apf::Numbering* numbers = apf::createNumbering(
       m,"my_numbers", m->getShape(), 1);
  //create an array to hold all our vertices
  apf::MeshEntity** vertices = new apf::MeshEntity*[total_elms];
  for(int counter = 0; counter < total_elms; counter++) {
    vertices[counter] = m->createVert(0);
  }
  //create a pointer to pass in quad vertices
  apf::MeshEntity* quad_verts[4];

  //rather than test for that zero,zero element every time
  //just instantiate it manually once
  apf::Vector3* temp_vec = new apf::Vector3(0,0,0);
  int vert_index = 0; // use this to place verts in right spot
  m->setPoint(vertices[vert_index],0,*(temp_vec));
  std::cout << vert_index << ": " << *(temp_vec) << std::endl;
  //now we can reuse by using Vector3::fromArray(const double* abc) 
  double xyz[3] = {0,0,0};
  for(int y_c = 1; y_c <= y_elms; y_c++) {
    //for the first column of x we make element to left
    xyz[0] = 0;
    xyz[1] = static_cast<double>(y_c) * y_size / 
      static_cast<double> (y_elms);
    temp_vec->fromArray(xyz); //load the values
    vert_index = (x_elms+1) * y_c;
    //print the location of each vertex
    std::cout << vert_index << ": " << *(temp_vec) << std::endl;
    m->setPoint(vertices[vert_index],0,*(temp_vec));
    for(int x_c = 1; x_c <= x_elms; x_c++) {
      xyz[0] = static_cast<double>(x_c) * x_size /
	  static_cast<double>(x_elms);
      //for the first row of y we construct an extra element below
      if(y_c == 1) {
	xyz[1] = 0;
	temp_vec->fromArray(xyz);
	vert_index = x_c;
	m->setPoint(vertices[vert_index],0,*(temp_vec));
	std::cout << vert_index << ": " << *(temp_vec) << std::endl;
      } //now create the actual row
      xyz[1] = static_cast<double>(y_c) * y_size / 
	static_cast<double>(y_elms);
      temp_vec->fromArray(xyz);
      vert_index = y_c*(x_elms+1) + x_c;
      m->setPoint(vertices[vert_index],0,*(temp_vec));
      //print the location of the vertex
      std::cout << vert_index << ": " << *(temp_vec) << std::endl;
      //create the quad element
      quad_verts[0] = vertices[vert_index];
      quad_verts[1] = vertices[(vert_index - 1)];
      quad_verts[2] = vertices[(vert_index - x_elms - 2)];
      quad_verts[3] = vertices[(vert_index - x_elms - 1)];
      //quad_verts[4] = vertices[vert_index];
      apf::Vector3 point;
      for(int p = 0; p < 4; p++) {
	m->getPoint(quad_verts[p], 0, point);
	std::cout << point  << std::endl;
      }
      apf::buildElement(m, 0, apf::Mesh::QUAD, quad_verts);
      std::cout << "quad created" << std::endl;
    }
  }
  //apf::deriveMdsModel(m);//this makes CAD model for classification
  //accept the changes
  //m->acceptChanges();
  //m->verify();

  //provided: after this should not change
  apf::writeVtkFiles("outQuad", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

