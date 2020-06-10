/***************************************************************************
 *   Copyright (C) 2006 by Jan Fostier   				   *
 *   jan.fostier@intec.ugent.be   					   *
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

#ifndef BMP_H
#define BMP_H

#include <cstdlib>
#include <iostream>

class Bmp {

public:
        /**
         * Load a bmp header from disc
         * @param width Width (in pixels) retreived from header (output)
         * @param height Height (in pixels) retreived from header (output)
         * @param filename File name of the input file
         */
        static void getBmpInfo(int *width, int *height,
                               const std::string &filename);

        /**
         * Write an uncompressed 24 bits bmp file to disc
         * @param data Pointer to a bitmap of data
         * @param resX Width of the bitmap (in pixels)
         * @param resY Height of the bitmap (in pixels)
         * @param filename File name of the output file
         */
        static void writeRGBbmp(void *data, int resX, int resY,
                                const std::string &filename);

        /**
         * Load an uncompressed 24 bits bmp file from disc
         * @param data Pointer to a bitmap of data
         * @param filename File name of the input file
         */
        static void loadRGBbmp(void *data, const std::string &filename);
};

#endif
