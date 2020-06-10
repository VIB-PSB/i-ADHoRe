/***************************************************************************
 *   Copyright (C) 2006 by Jan Fostier                                     *
 *   jan.fostier@intec.ugent.be                                            *
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

#ifndef GRAFIX_H
#define GRAFIX_H

#include <cstdlib>
#include <iostream>

#include "color.h"
#include "pnghandling.h"

// ========================================================================
// GRAFIX
// ========================================================================

class Grafix {

public:
        /**
         * Create a resX * resY canvas
         * @param resX x-resolution
         * @param resY y-resolution
         */
        Grafix(int resX, int resY);

        /**
         * Load a bmp file from disc
         * @param filename filename of the bmp file
         */
        Grafix(const std::string &filename);

        /**
         * Destructor
         */
        ~Grafix();

        /**
         * Get the resolution
         * @param resX x-resolution (output)
         * @param resY y-resolution (output)
         */
        void getResolution(int &resX, int &resY);

        /**
         * Set the current drawing color
         * @param c
         */
        void setDrawingColor(const Color &c);

        /**
         * Calculate viewport
         * @param LLx Lower-left x-coordinate
         * @param LLy Lower-left y-coordinate
         * @param URx Upper-right x-coordinate
         * @param URy Upper-right y-coordinate
         */
        void setCanvasViewport(double LLx, double LLy, double URx, double URy);

        /**
         * Save the canvas to a bmp file
         * @param filename
         */
        void saveCanvasBmp(const std::string &filename);

        #ifdef HAVE_PNG
        void saveCanvasPng(const std::string &filename);
        #endif
        // ===================================================================
        // PRIMITIVE DRAWING FUNCTIONS (NO VIEWPORT TRANSFORMATION)
        // ===================================================================

        /**
         * Put a pixel
         * @param x x-coordinate
         * @param y y-coordinate
         */
        void putPixel(int x, int y);

        /**
         * Get the color in the bitmap
         * @param x x-coordinate
         * @param y x-coordinate
         * @return Color
         */
        Color getPixel(int x, int y);

        /**
         * Draw a line using the Bresenham algorithm
         * @param x1 x-coordinate p1
         * @param y1 y-coordinate p1
         * @param x2 x-coordinate p2
         * @param y2 y-coordinate p2
         */
        void lineBresenham(int x1, int y1, int x2, int y2);

        void drawBox(int begin_x, int end_x, int begin_y, int end_y);

        /**
         * Draw a line using the Bresenham algorithm
         * @param x1 x-coordinate p1
         * @param y1 y-coordinate p1
         * @param x2 x-coordinate p2
         * @param y2 y-coordinate p2
         */
        void floatLine(float x1, float y1, float x2, float y2);

        /**
         * Draw a circle using the Bresenham algorithm
         * @param x x-coordinate centre
         * @param y y-coordinate centre
         * @param radius radius of the circle
         */
        void circleBresenham(int x, int y, int radius);

        // ===================================================================
        // PRIMITIVE DRAWING FUNCTIONS (WITH VIEWPORT TRANSFORMATION)
        // ===================================================================

        /**
         * Draw a pixel on the bitmap (with viewport TF)
         * @param x x-coordinate
         * @param y y-coordinate
         */
        void drawPixel(double x, double y);

        /**
         * Draw a line using the Bresenham algorithm (with viewport TF)
         * @param x1 x-coordinate p1
         * @param y1 y-coordinate p1
         * @param x2 x-coordinate p2
         * @param y2 y-coordinate p2
         */
        void drawLine(double x1, double y1, double x2, double y2);

        /**
         * Draw a circle using the Bresenham algorithm (with viewport TF)
         * @param x x-coordinate centre
         * @param y y-coordinate centre
         * @param radius radius of the circle
         */
        void drawCircle(double x, double y, double radius);

        /**
         * Floodfill the bitmap
         * @param x x-coordinate startpoint
         * @param y y-coordinate startpoint
         */
        void seedFill(double x, double y);

        // ===================================================================
        // AUXILIARY FUNCTIONS
        // ===================================================================

        /**
         * Calculate viewport transformation of a world coordinate
         * @param x x-coordinate (world)
         * @param y y-coordinate (world)
         * @param x_o x-coordinate (screen - output)
         * @param y_o y-coordinate (screen - output)
         */
        void viewportTF(double x, double y, int &x_o, int &y_o);

        /**
         * Calculate viewport transformation of a screen coordinate
         * @param x x-coordinate (screen)
         * @param y y-coordinate (screen)
         * @param x_o x-coordinate (world - output)
         * @param y_o y-coordinate (world - output)
         */
        void inverseViewportTF(int x, int y, double &x_o, double &y_o);

private:

        /**
         * Plots all eight quadrants of a circle
         * @param x x-coordinate of the centre
         * @param y y-coordinate of the centre
         * @param x1 x-coordinate (screen) to plot
         * @param y1 y-coordinate (screen) to plot
         */
        void plotCircle(int x, int y, int x1, int y1);

        /**
         * Recursive flood fill
         * @param x x-coordinate (screen)
         * @param y y-coordinate (screen)
         */
        void seedFillRecursive(int x, int y);

        int resX;               // canvas dimension X
        int resY;               // canvas dimension Y
        Color *canvas;          // canvas pointer
        Color drawCol;          // drawing color
        Color ffBg;             // floodfill background color

        // viewport transformation
        double vp_scale;
        double vp_off_x;
        double vp_off_y;
};

#endif
