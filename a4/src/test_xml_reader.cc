#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/util/XMLString.hpp>

#include "vtkSAX2Handler.h"

#include <iostream>
#include <cstdio>

using namespace std;

int main (int argc, char* argv[]) {

    if(argc != 2) {
        printf("Usage %s filename\n",argv[0]);
        exit(0);
    }
    try {
        XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Error during initialization! :\n";
        cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release(&message);
        return 1;
    }

    std::cout << "file target is " << argv[1] << std::endl;

    SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, true);   // optional

    //DefaultHandler* defaultHandler = new DefaultHandler();
    vtkSAX2Handler* defaultHandler = new vtkSAX2Handler();
    parser->setContentHandler(defaultHandler);
    parser->setErrorHandler(defaultHandler);

    try {
        parser->parse(argv[1]);
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release(&message);
        return -1;
    }
    catch (const SAXParseException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release(&message);
        return -1;
    }
    catch (...) {
        cout << "Unexpected Exception \n" ;
        return -1;
    }

    delete parser;
    delete defaultHandler;
    /*clear all mememory useage*/
    XMLPlatformUtils::Terminate();
    return 0;
}
