#include "SvgWriter.h"

typedef unsigned int uint;

SvgWriter::SvgWriter(const std::string filename, double xDim, double yDim) : xsize(xDim), ysize(yDim)
{
    ofs.open(filename.c_str());
    writeHeader();
}

SvgWriter::~SvgWriter()
{
    if (ofs.is_open())
        finish();
}

void SvgWriter::finish()
{
    ofs << "</svg>";
    ofs << endl;
    ofs.close();
}

string SvgWriter::convertColor(const Color& color)
{
    char buffer[8];
    uint r,g,b;
    color.getIntegerRGB(r,g,b);

    sprintf(buffer,"#%x%x%x",(unsigned char)r, (unsigned char)g, (unsigned char)b);
    string colorToString(buffer);
    return colorToString;
}

void SvgWriter::drawRectangle(double x, double y, double width, double height, const Color& color)
{
    string fill=convertColor(color);
    string stroke="black";
    ofs << "<rect ";
    setXMLLengthVariable(ofs,x)
    setXMLLengthVariable(ofs,y)
    setXMLLengthVariable(ofs,width)
    setXMLLengthVariable(ofs,height)
    setXMLVariable(ofs,fill)
    //setXMLVariable(ofs,stroke)
    ofs << "stroke-width=\"" << 5 << "\" ";
    closeXMLVariable(ofs)
}



void SvgWriter::drawDot(double cx, double cy, double r)
{
    ofs << "<circle ";
    setXMLLengthVariable(ofs,cx)
    setXMLLengthVariable(ofs,cy)
    setXMLLengthVariable(ofs,r)
    closeXMLVariable(ofs)

}

void SvgWriter::drawLine(double x1, double y1, double x2, double y2)
{
    string stroke="black";
    ofs << "<line ";
    setXMLLengthVariable(ofs,x1)
    setXMLLengthVariable(ofs,y1)
    setXMLLengthVariable(ofs,x2)
    setXMLLengthVariable(ofs,y2)
    setXMLVariable(ofs,stroke)
    //ofs << "stroke=\"" << stroke << "\" ";
    ofs << "stroke-width=\"" << 2 << "\" ";
    closeXMLVariable(ofs)
}

void SvgWriter::drawText(double x, double y, string textString, int fontsize)
{
    string fill="black";
    ofs << "<text ";
    setXMLLengthVariable(ofs,x)
    setXMLLengthVariable(ofs,y)
    ofs << "font-family=\"Verdana\" " <<"font-size=\"" << fontsize << "\" ";
    setXMLVariable(ofs,fill)
    ofs << ">" << endl;
    ofs << textString << endl;
    ofs << "</text>"  << endl;
}

void SvgWriter::writeHeader()
{

    double width=xsize;
    double height=ysize;
    ofs << "<?xml version=\"1.0\" standalone=\"yes\"?>";
    ofs << endl;

    ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
    setXMLLengthVariable(ofs,width)
    setXMLLengthVariable(ofs,height)
    ofs << ">" << endl;
}
