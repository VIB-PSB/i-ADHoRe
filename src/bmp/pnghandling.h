#ifndef PNGHANDLING_H
#define PNGHANDLING_H

#ifdef HAVE_PNG

#include <png.h>
#include "color.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <string>
using std::string;


//NOTE: example code used and modified from http://zarb.org/~gc/html/libpng.html
/* Copyright 2002-2008 Guillaume Cottenceau.
* This software may be freely redistributed under the terms
* of the X11 license.
*/

//NOTE For advanced use extra parameters can be set, this class was written to allow for elementary uses
//Information to expand this class can be found in libpng.txt written by the PNG Development Group

class PngHandling
{
public:

    /**
    * Constructor
    */
    PngHandling();

    /**
    * Destructor
    */
    ~PngHandling();

    /**
    * @return true if png_ptr and info_ptr are succesfully constructed
    */
    bool isInitialized() const;

    /**
    * @param w number of pixels in x-direction
    * @param h number of pixels in y-direction
    */
    void setImageDimensions(uint w, uint h);

    /**
    * Allocate memory for PNG image & initialize background
    */
    void prepareImage();

    void putPixel(uint x, uint y, const Color& c);

    void writeToFile(char* file_name);

protected:

    int width, height;
    png_byte color_type, bit_depth;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep * row_pointers;

    void putBackground(const Color& c);

};

#endif
#endif
