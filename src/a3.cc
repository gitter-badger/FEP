#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <stdio.h>
#include <stdlib.h>
#include <fstream>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("cube.dmg", "parallelMesh.smb");
  
  //test the idea of multi process
  if(PCU_Comm_Self()==0) {
    std::cout << "there are " << PCU_Comm_Peers() << "total processes"
	      << std::endl;    
  }
  std::cout << "hello from process " << PCU_Comm_Self() << std::endl;
  int i = PCU_Comm_Self() + 1;
  PCU_Add_Ints(&i, 1);
  int b[2] = {i+1, i+4};
  PCU_Add_Ints(b,2);
  std::cout << i << b[0] << b[1] << std::endl;
  

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
