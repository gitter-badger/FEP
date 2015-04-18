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

  char* message;
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  apf::Adjacent up;

  if(PCU_Comm_Self() == 0) {
    std::cout << "number of peers" << PCU_Comm_Peers() << std::endl;
    it = m->begin(0);
    apf::Copies remotes;
    while(e = m->iterate(it)) {
      m->getRemotes(e, remotes);
      for(apf::Copies::iterator it = remotes.begin(); it != remotes.end(); ++it) {
	std::cout << "shared with:" << it->second << std::endl; 
      }
    }
  }
  
  apf::Numbering* vert_nums = apf::numberOwnedDimension(m, "verts", 0);
  
 
  apf::writeVtkFiles("part",m); 
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
