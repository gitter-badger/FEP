#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

/*store the mesh pointer we are building directly*/
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfNumberingClass.h>
#include <apfShape.h>
#include <apfFieldData.h>

using namespace xercesc_3_1;

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

    void characters(
        const   XMLCh* const    chars,
        const   XMLSize_t       length
    );

    void error(const SAXParseException &exc);

    void fatalError(const SAXParseException&);

    apf::Mesh* mesh;
private:
    bool _have_seen_VTK_file_tag;
    bool _ignore_tag;
    bool _ignore_tag_char_data;
};
