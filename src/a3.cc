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

  apf::Numbering* vert_nums = apf::numberOwnedDimension(m, "verts", 0);
  apf::Numbering* reg_nums = apf::numberOwnedDimension(m, "regions", 3);

  if(PCU_Comm_Self() == 0) {
    std::cout << "number of peers" << PCU_Comm_Peers() << std::endl;
    it = m->begin(0);
    apf::Copies remotes;
    while(e = m->iterate(it)) {
      m->getRemotes(e, remotes);
      if(m->isOwned(e)) {
	for(apf::Copies::iterator it = remotes.begin(); it != remotes.end(); ++it) {
	  std::cout << e << " shared with:" << it->second << std::endl; 
	}
      }
      std::cout << "====loop====" << std::endl;
    }
  }
  //Part 3
  apf::writeVtkFiles("pre_migr", m);
  //migration plan should exist for each process
  apf::Migration* plan = new apf::Migration(m);
  
  //get the mesh regions that are on boundary
  it = m->begin(3);
  //count the number of faces on the partition model face
  int partition_face_cnt = 0;
  while(e = m->iterate(it)) {
    //we assume that there are no shadow regions
    //otherwise the adjacencies would not make much sense
    apf::Downward faces;
    int num_faces = m->getDownward(e, 2, faces);
    int migrate_flag = 0;
    int destination = -1;
    for(int ii = 0; ii < num_faces; ++ii) {
      if( ! m->isOwned(faces[ii])) {
	++partition_face_cnt;
	if(PCU_Comm_Self() == 0) { 
	  if(! plan->has(e)) {
	    plan->send(e, m->getOwner(faces[ii]));
	  }
	}
      }
    }
  }

  printf("Partition %d reports %d partition faces before migration\n", \
	 PCU_Comm_Self() , partition_face_cnt);
  //now send all of the migration plan
  m->migrate(plan);

  //now compute the number of shared faces
  it = m->begin(3);
  //count the number of faces on the partition model face
  int aftr_partition_face_cnt = 0;
  while(e = m->iterate(it)) {
    //we assume that there are no shadow regions
    //otherwise the adjacencies would not make much sense
    apf::Downward faces;
    int num_faces = m->getDownward(e, 2, faces);
    int migrate_flag = 0;
    int destination = -1;
    for(int ii = 0; ii < num_faces; ++ii) {
      if( ! m->isOwned(faces[ii])) {
	++aftr_partition_face_cnt;
      }
    }
  }
  printf("Partition %d reports %d partition faces after migration\n", \
	 PCU_Comm_Self() , aftr_partition_face_cnt);

 
  apf::writeVtkFiles("post_migr", m); 
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
