#include "vtkSAX2Handler.h"
#include <iostream>


using namespace std;
using namespace xercesc_3_1;

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
}

void vtkSAX2Handler::fatalError(const SAXParseException& exception)
{
    char* message = XMLString::transcode(exception.getMessage());
    cout << "Fatal Error: " << message
         << " at line: " << exception.getLineNumber()
         << endl;
    XMLString::release(&message);
}
