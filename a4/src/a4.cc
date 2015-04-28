#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <gtest/gtest.h>

int main(int argc, char** argv)
{
	MPI_Init(&argc,&argv);
  	PCU_Comm_Init();
  	::testing::InitGoogleTest(&argc, argv);
	int rc = RUN_ALL_TESTS();
	if(rc != 0){
		printf("Google Test returned: %d", rc);
	}
	PCU_Comm_Free();
  	MPI_Finalize();
	return 0;
}
