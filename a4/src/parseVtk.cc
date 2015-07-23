#include <fstream>
#include <sstream>

#include <xercesc/sax2/DefaultHandler.hpp>

#include <apf.h>
/*apf.h includes the definition of the writeVtkFiles,
note that most of the helper functions called by
writeVtkFiles are static, which is no good to us,
maybe we will ask for some reuse? */
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfNumberingClass.h>
#include <apfShape.h>
#include <apfFieldData.h>


static void readPvtuFile(const char* filename, apf::Mesh* m) {
	//apf::writeVtkFiles("foo", m);
}

static void readVtuFile(const char* filename, apf::Numbering* n) {
	

	
}
