#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

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
   		const 	XMLCh* const uri,
   		const 	XMLCh* const localname,
   		const 	XMLCh* const qname
   	);

    void fatalError(const SAXParseException&);
};
