/***************************************************************************
 *   Copyright (C) 2006 by Jan Fostier   				   *
 *   jan.fostier@intec.ugent.be   			           	   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "bmp.h"
#include <fstream>

using namespace std;

// necessary to force 2 byte aligning en lieu de 4 or 8
#pragma pack(2)

// ========================================================================
// BMP FILES
// ========================================================================

// format specific structures
typedef struct {
        unsigned short int type;
        unsigned int size;
        unsigned short int reserved1, reserved2;
        unsigned int offset;
} HEADER;

typedef struct {
        unsigned int size;
        int width,height;
        unsigned short int planes;
        unsigned short int bits;
        unsigned int compression;
        unsigned int imagesize;
        int xresolution,yresolution;
        unsigned int ncolours;
        unsigned int importantcolours;
} INFOHEADER;

void Bmp::writeRGBbmp(void *data, int resX, int resY,
                      const std::string &filename)
{
        INFOHEADER infoheader;
        HEADER header;

        header.type = 19778;            // water brand for bmp files
        header.size = sizeof(header)+sizeof(infoheader)+resX*resY*3;
        header.reserved1 = 0;
        header.reserved2 = 0;
        header.offset = 54;

        infoheader.size = sizeof(infoheader);
        infoheader.width = resX;
        infoheader.height = resY;
        infoheader.planes = 1;
        infoheader.bits = 24;           // 24 bits per pixel
        infoheader.compression = 0;     // no compression
        infoheader.imagesize = 0;       // allowed 0 because no compression
        infoheader.xresolution = 0;
        infoheader.yresolution = 0;
        infoheader.ncolours = 0;
        infoheader.importantcolours = 0;

        // open file and write data and header
        ofstream ofs;
        ofs.open(filename.c_str(), ios::binary);

        ofs.write(reinterpret_cast<char*>(&header), sizeof(header));
        ofs.write(reinterpret_cast<char*>(&infoheader), sizeof(infoheader));
        ofs.write(reinterpret_cast<char*>(data), 3*resX*resY);
        ofs.close();
}

void Bmp::getBmpInfo(int *width, int *height, const std::string &filename)
{
        INFOHEADER infoheader;
        HEADER header;

        // open file and load header
        ifstream ifs;
        ifs.open(filename.c_str(), ios::binary);

        ifs.read(reinterpret_cast<char*>(&header), sizeof(header));
        ifs.read(reinterpret_cast<char*>(&infoheader), sizeof(infoheader));

        ifs.close();

        *width = infoheader.width;
        *height = infoheader.height;
}

void Bmp::loadRGBbmp(void *data, const std::string &filename)
{
        INFOHEADER infoheader;
        HEADER header;

        // open file and load header
        ifstream ifs;
        ifs.open(filename.c_str(), ios::binary);

        ifs.read(reinterpret_cast<char*>(&header), sizeof(header));
        ifs.read(reinterpret_cast<char*>(&infoheader), sizeof(infoheader));
        ifs.read(reinterpret_cast<char*>(data),
                 3*infoheader.width*infoheader.height);

        ifs.close();
}
