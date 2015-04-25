#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <gtest/gtest.h>

int main(int argc, char const *argv[])
{
	printf("Begin the studpid\n");
	int rc = RUN_ALL_TESTS();
	printf("%d\n",rc );
	return 0;
}