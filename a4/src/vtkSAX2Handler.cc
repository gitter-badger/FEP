#include "vtkSAX2Handler.h"
#include <iostream>
#include <cstring>

#define VTKFILE_TAG_NAME  "VTKFile"

using namespace std;

vtkSAX2Handler::vtkSAX2Handler()
{
    this->_have_seen_VTK_file_tag = false;
    this->_ignore_tag = false;
}

void vtkSAX2Handler::startElement(const   XMLCh* const    uri,
                            const   XMLCh* const    localname,
                            const   XMLCh* const    qname,
                            const   Attributes&     attrs)
{
    char* message = XMLString::transcode(localname);
    cout << "I saw element: "<< message << endl;

    /*check for VTKFile tag and start */
    if( std::strcmp(message, VTKFILE_TAG_NAME ) == 0) {
        if(this->_have_seen_VTK_file_tag == true) {
            std::cout << "\tSEEN TWICE" << std::endl;
            /*two nested VTKFile tags is an error, should never happen*/
            XMLCh *error_message, *publicId, *systemId;
            XMLFileLoc lineNumber, columnNumber;
            const char* error_msg = "Nested <VTKFile> tag not allowed";
            error_message = XMLString::transcode(error_msg);
            publicId = XMLString::transcode("a publicId");
            systemId = XMLString::transcode("a systemID");
            /*right now we don't know where in the document we are*/
            lineNumber = 0;
            columnNumber = 0;
            /*leaving off mememory manager*/
            throw SAXParseException(error_message, publicId, systemId, lineNumber, columnNumber);
        }
        this->_have_seen_VTK_file_tag = true;
        XMLString::release(&message);
    } else {
        /*if the tage is not a VTKFile then do not bother parsing rest*/
        XMLString::release(&message);
        return;
    }

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

    if(std::strcmp(message, VTKFILE_TAG_NAME) == 0) {
        this->_have_seen_VTK_file_tag = false;
    }

    XMLString::release(&message);
}

void vtkSAX2Handler::characters(
        const   XMLCh* const    chars,
        const   XMLSize_t       length
    ) {

    std::cout << "character data of length " << length << std::endl;
}

void vtkSAX2Handler::fatalError(const SAXParseException& exception)
{
    char* message = XMLString::transcode(exception.getMessage());
    cout << "Fatal Error: " << message
         << " at line: " << exception.getLineNumber()
         << endl;
    XMLString::release(&message);
}

