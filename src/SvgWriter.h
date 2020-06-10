#ifndef SVGWRITER_H
#define SVGWRITER_H

#include <iostream>
#include <fstream>
#include "bmp/color.h"

using namespace std;

#define setXMLVariable(streamer,x)  streamer << #x <<  "=\"" << x << "\" ";
#define setXMLLengthVariable(streamer,x) streamer << #x <<  "=\"" << x << "cm\" ";
#define closeXMLVariable(streamer)  streamer << "/>" << endl;


//Resource1: http://www.w3.org/TR/SVG

//Resource2: Getting Started with XML: A Manual and Workshop
// by Eric Lease Morgan


class SvgWriter
{

public:
    /**
    * Constructor
    @param filename name of the xml file
    */
    SvgWriter(const std::string filename, double xDim, double yDim);

    /**
    * Destructor
    */
    ~SvgWriter();

    /**
    * Closes XML elements and writes file
    */
    void finish();

    /**
    * Convert color to a string which can be used by the XML macros
    */
    string convertColor(const Color& color);

    /**
    * Draws a rectangle
    *@param x x-coordinate of upper left corner point
    *@param y y-coordinate of upper left corner point
    */
    void drawRectangle(double x, double y, double width, double height, const Color& color);

    /**
    *@param cx center of circle x-coordinate
    *@param cy center of circle y-coordinate
    *@param r radius of dot
    */
    void drawDot(double cx, double cy, double r);

    void drawLine(double x1, double y1, double x2, double y2);

    void drawText(double x, double y, string textString, int fontsize);

private:

    double xsize,ysize;
    ofstream ofs;

    SvgWriter(){};
    SvgWriter(const SvgWriter& sw){};
    void operator=(const SvgWriter& sw){};

    /**
    * Write XML header
    */
    void writeHeader();

};

#endif
