#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

/*store the mesh pointer we are building directly*/
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfNumberingClass.h>
#include <apfShape.h>
#include <apfFieldData.h>


class vtkSAX2Handler : public DefaultHandler {
public:
    vtkSAX2Handler(); /*public constructor*/

    void startElement(
        const   XMLCh* const    uri,
        const   XMLCh* const    localname,
        const   XMLCh* const    qname,
        const   Attributes&     attrs
    );

    void endElement(
        const   XMLCh* const uri,
        const   XMLCh* const localname,
        const   XMLCh* const qname
    );

    void fatalError(const SAXParseException&);

    apf::Mesh* mesh;

};
