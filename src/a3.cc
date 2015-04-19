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
#include <vector>

#define INT_TAG_LENGTH 3


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh("cube.dmg", "parallelMesh.smb");

  //Part 1 and 2
  //iterate through each vertex and compute adjacency for each
  apf::Numbering* vert_nums = apf::numberOwnedDimension(m, "verts", 0);
  apf::Numbering* edge_nums = apf::numberOwnedDimension(m, "edges", 1);
  apf::Numbering* face_nums = apf::numberOwnedDimension(m, "faces", 2);
  apf::Numbering* reg_nums = apf::numberOwnedDimension(m, "regions", 3);

  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* vert;
  apf::MeshEntity* e;
  apf::Adjacent up;
  apf::MeshTag* local_adj_tag = m->createIntTag("all_adj", 3);
  //iterate over each mesh vertex and get all adjacencies
  int current_partition = PCU_Comm_Self();
  int efr_adj[INT_TAG_LENGTH];
  //use tag size to record  with number of bytes in our tag
  int TRANSFER_TAG_SIZE = sizeof(int) *INT_TAG_LENGTH;

  PCU_Comm_Begin();
  std::cout << "here" << std::endl;
  while(vert = m->iterate(it)) {
    std::cout << "vert: " << vert << " " << current_partition <<  std::endl;
    for(int dim = 1; dim < 4; ++dim) {
      std::cout << "dimension " << dim << current_partition << std::endl;
      m->getAdjacent(vert, dim, up);
      std::cout << "got adjacent " << current_partition << std::endl; 
      int owned_dim_count = 0;
      for(int ii = 0; ii < up.getSize(); ++ii) {
	//only count adjacencies on the local part
	//this avoid double counting
	if(m->isOwned(up[ii])) {
	  ++owned_dim_count;
	}
      }
      //offset by one
      int index = dim - 1;
      efr_adj[index] = owned_dim_count;
    }
    std::cout << "===loop" << current_partition << "===" << std::endl;
    //now add the data to the tag
    m->setIntTag(vert, local_adj_tag, efr_adj);
    std::cout << "tag set " << current_partition << "bool " << m->isOwned(vert) <<  std::endl;
    //now send the taged data to the other partitions
    //only send the data that other partitions want
    if( ! m->isOwned(vert)) {
      std::cout << "not owned " << current_partition << std::endl;
      //send the vertex first
      int target_process = m->getOwner(vert);
      std::cout << sizeof(apf::MeshEntity*) << " " << target_process  << current_partition << std::endl;
      int rc = PCU_Comm_Write(target_process, (void*)(&vert), sizeof(apf::MeshEntity*));
      std::cout << "sent pointer" << current_partition << std::endl;
      PCU_Comm_Write(target_process, reinterpret_cast<void*>(&efr_adj), TRANSFER_TAG_SIZE );
      std::cout << "sent data" << current_partition << std::endl;
    }
  }
  std::cout << "boop from: " << current_partition << std::endl; 
  //we are done with local adjacencies
  PCU_Comm_Send();
  std::cout << "done sending" << std::endl;
  while(PCU_Comm_Receive()) {
    void* data;
    PCU_Comm_Unpack(data, sizeof(apf::MeshEntity*));
    apf::MeshEntity* target_vert =  reinterpret_cast<apf::MeshEntity*>(data);
    //parts should send adjacency to correct part number
    assert(m->isOwned(target_vert));

    PCU_Comm_Unpack(data,TRANSFER_TAG_SIZE );
    int foreign_adj[INT_TAG_LENGTH];
    memcpy(reinterpret_cast<void*>(foreign_adj), data, TRANSFER_TAG_SIZE );
    //get the current adj tag and add
    int temp_tag[INT_TAG_LENGTH];
    m->getIntTag(target_vert, local_adj_tag, temp_tag);
    //sum the current adjacencies and the local adjacencies
    for(int jj = 0; jj < INT_TAG_LENGTH; ++jj) {
      temp_tag[jj] += foreign_adj[jj];
    } 
    m->setIntTag(target_vert, local_adj_tag, temp_tag);
  }
  
  //now report the adjacencies for each vertex in the local part
  it = m->begin(0);
  while(vert = m->iterate(it)) {
    if(m->isOwned(vert)) {
      int temp_tag[INT_TAG_LENGTH];
      m->getIntTag(vert, local_adj_tag, temp_tag);
      printf("L: %d ZVert num: %d, E: %d F: %d R: %d", current_partition, \
	     apf::getNumber(vert_nums, vert, 0, 0), temp_tag[0], temp_tag[1], temp_tag[2]);
    }
  }
  
  

  //Part 3
  apf::writeVtkFiles("pre_migr", m);
  //migration plan should exist for each process
  apf::Migration* plan = new apf::Migration(m);
  
  //get the mesh regions that are on boundary
  it = m->begin(2);
  apf::Adjacent up_r;
  //count the number of faces on the partition model face
  int partition_face_cnt = 0;
  while(e = m->iterate(it)) {
    //we assume that there are no shadow regions
    //otherwise the adjacencies would not make much sense
    if( ! m->isOwned(e)) {
      m->getAdjacent(e, 3, up_r); 
      assert(up_r.getSize() == 1);
      ++partition_face_cnt;
      if(PCU_Comm_Self() == 0) { 
	apf::MeshEntity* region = up_r[0];
	if(! plan->has(region)) {
	  plan->send(region, m->getOwner(e));
	}
       }
    }
  }

  printf("Partition %d reports %d partition faces before migration\n", \
	 PCU_Comm_Self() , partition_face_cnt);
  //now send all of the migration plan
  m->migrate(plan);

  //now compute the number of shared faces
  it = m->begin(2);
  //count the number of faces on the partition model face
  int aftr_partition_face_cnt = 0;
  while(e = m->iterate(it)) {
    //we assume that there are no shadow regions
    //otherwise the adjacencies would not make much sense
    if( ! m->isOwned(e)) {
      m->getAdjacent(e, 3, up_r); 
      assert(up_r.getSize() == 1);
      ++aftr_partition_face_cnt;
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



