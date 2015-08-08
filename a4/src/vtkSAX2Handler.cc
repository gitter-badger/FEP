#include "vtkSAX2Handler.h"
#include <iostream>
#include <cstring>
#include <string>

#define VTKFILE_TAG_NAME  "VTKFile"
#define DATA_ARRAY_TAG_NAME "DataArray"

using namespace std;

vtkSAX2Handler::vtkSAX2Handler()
{
    this->_have_seen_VTK_file_tag = false;
    this->_ignore_tag = false;
    /*only one tag will have character data, so we turn
    * parsing of any character data off as a sensible default*/
    this->_ignore_tag_char_data = true;
}

void vtkSAX2Handler::startElement(const   XMLCh* const    uri,
                            const   XMLCh* const    localname,
                            const   XMLCh* const    qname,
                            const   Attributes&     attrs)
{
    char* message = XMLString::transcode(localname);
    std::string tag_name(message);
    XMLString::release(&message);

    cout << "I saw element: "<< tag_name << endl;

    /*check for VTKFile tag and start */
    if(tag_name == std::string(VTKFILE_TAG_NAME)) {
        this->_have_seen_VTK_file_tag = true;
    }
    if(this->_have_seen_VTK_file_tag == false) {
        /*no need to process anything */
        return;
    }

    /*only accept character data if we are in a data array tag*/
    if(tag_name == std::string(DATA_ARRAY_TAG_NAME)) {
        this->_ignore_tag_char_data = false;
    }


    // XMLSize_t len, ii;
    // len = attrs.getLength();
    // std::cout << "  number of attrs: = " << len << std::endl;
    // const XMLCh* attr_name;
    // for(ii = 0; ii < len; ++ii){
    //     attr_name = attrs.getQName(ii);
    //     message = XMLString::transcode(attr_name);
    //     std::cout << "attr: " << message << " value: " \
    //         << XMLString::transcode(attrs.getValue(ii)) \
    //         << std::endl;
    //     XMLString::release(&message); 
    // }
}

void vtkSAX2Handler::endElement(const   XMLCh* const uri,
        const   XMLCh* const localname,
        const   XMLCh* const qname)
{
    char* message = XMLString::transcode(localname);
    std::string tag_name(message);
    XMLString::release(&message);
    std::cout << "\tend of element: " << tag_name << std::endl;

    if(tag_name == std::string(VTKFILE_TAG_NAME)) {
        this->_have_seen_VTK_file_tag = false;
    }

    /*only accept character data if we are in a data array tag*/
    if(tag_name == std::string(DATA_ARRAY_TAG_NAME)) {
        this->_ignore_tag_char_data = false;
    }

}

void vtkSAX2Handler::characters(
        const   XMLCh* const    chars,
        const   XMLSize_t       length
    ) {
    if(this->_ignore_tag_char_data) {
        return;
    }
    std::cout << "character data of length " << length << std::endl;
}

void vtkSAX2Handler::error(const SAXParseException& exception) {
    /*right now just throw on all errors*
    * all of the DTD errors use this callback, and right
    * now I have no way to distinguish them, so throw on all
    */    
    throw SAXParseException(exception);
}

void vtkSAX2Handler::fatalError(const SAXParseException& exception)
{
    char* message = XMLString::transcode(exception.getMessage());
    cout << "Fatal Error: " << message
         << " at line: " << exception.getLineNumber()
         << endl;
    XMLString::release(&message);
}

