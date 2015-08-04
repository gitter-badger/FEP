#include "vtkSAX2Handler.h"
#include <iostream>


using namespace std;

vtkSAX2Handler::vtkSAX2Handler()
{
}

void vtkSAX2Handler::startElement(const   XMLCh* const    uri,
                            const   XMLCh* const    localname,
                            const   XMLCh* const    qname,
                            const   Attributes&     attrs)
{
    char* message = XMLString::transcode(localname);
    cout << "I saw element: "<< message << endl;
    XMLString::release(&message);
    XMLSize_t len, ii;
    len = attrs.getLength();
    std::cout << "  number of attrs: = " << len << std::endl;
    const XMLCh* attr_name;
    for(ii = 0; ii < len; ++ii){
        attr_name = attrs.getQName(ii);
        message = XMLString::transcode(attr_name);
        std::cout << "attr: " << message << " value: " \
            << XMLString::transcode(attrs.getValue(ii)) \
            << std::endl;
        XMLString::release(&message); 
    }
}

void vtkSAX2Handler::endElement(const   XMLCh* const uri,
        const   XMLCh* const localname,
        const   XMLCh* const qname)
{
    char* message = XMLString::transcode(localname);
    std::cout << "\tend of element: " << message << std::endl;
    XMLString::release(&message);
}

void vtkSAX2Handler::fatalError(const SAXParseException& exception)
{
    char* message = XMLString::transcode(exception.getMessage());
    cout << "Fatal Error: " << message
         << " at line: " << exception.getLineNumber()
         << endl;
    XMLString::release(&message);
}

