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

#include "grafix.h"
#include "bmp.h"
#include <cassert>
#include <math.h>

// ========================================================================
// GRAFIX
// ========================================================================

// create a canvas with specified dimensions
Grafix::Grafix(int resX_, int resY_) : resX(resX_), resY(resY_), vp_scale(1),
               vp_off_x(0), vp_off_y(0)
{
        assert(resX > 0);
        assert(resY > 0);
        assert(resX%4==0);

        canvas = new Color[resX*resY];
}

Grafix::Grafix(const std::string &filename) : vp_scale(1), vp_off_x(0),
               vp_off_y(0)
{
        Bmp::getBmpInfo(&resX, &resY, filename);
        canvas = new Color[resX*resY];
        Bmp::loadRGBbmp(canvas, filename);
}

Grafix::~Grafix()
{
        delete [] canvas;
}

void Grafix::getResolution(int &resX_, int &resY_)
{
        resX_ = resX;
        resY_ = resY;
}

void Grafix::setDrawingColor(const Color &c)
{
        drawCol = c;
}

void Grafix::setCanvasViewport(double LLx, double LLy, double URx, double URy)
{
        double dx = URx - LLx;
        double dy = URy - LLy;

        double ax = (double)resX/dx;
        double ay = (double)resY/dy;

        if (ax < ay) {
                vp_scale = ax;
                vp_off_x = -vp_scale*LLx;
                vp_off_y = (double)resY/2.0-vp_scale*(LLy+dy/2.0);
        } else {
                vp_scale = ay;
                vp_off_x = (double)resX/2.0-vp_scale*(LLx+dx/2.0);
                vp_off_y = -vp_scale*LLy;
        }
}

void Grafix::saveCanvasBmp(const std::string &filename)
{
        Bmp::writeRGBbmp(canvas, resX, resY, filename);
}

#ifdef HAVE_PNG
void Grafix::saveCanvasPng(const std::string &filename)
{
    PngHandling png;
    png.setImageDimensions(resX, resY);
    png.prepareImage();
    for (int x=0; x<resX; x++) {
        for (int y=0; y<resY; y++) {
            png.putPixel(x,y,getPixel(x,y));
        }
    }
    png.writeToFile(const_cast<char*>(filename.c_str()));
}
#endif

// ========================================================================
// PRIMITIVE DRAWING FUNCTIONS (NO VIEWPORT TRANSFORMATION)
// ========================================================================

void Grafix::putPixel(int x, int y)
{
        if ((x >= 0) && (x < resX) && (y >= 0) && (y < resY))
                canvas[y*resX+x] = drawCol;
}

Color Grafix::getPixel(int x, int y)
{
        return canvas[y*resX+x];
}

void Grafix::lineBresenham(int x1, int y1, int x2, int y2)
{
        // offset and stepping
        int dy = y2-y1;
        int dx = x2-x1;
        int stepx, stepy;

        if (dy < 0) { dy = -dy;  stepy = -1; } else { stepy = 1; }
        if (dx < 0) { dx = -dx;  stepx = -1; } else { stepx = 1; }
        dy <<= 1; // dy is now 2*dy
        dx <<= 1; // dx is now 2*dx

        // put first pixel
        putPixel(x1, y1);

        if (dx > dy) { // slope is less than 45 degrees
                int fraction = dy - (dx >> 1);  // same as 2*dy - dx
                while (x1 != x2) {
                        if (fraction >= 0) {
                                y1 += stepy;
                                fraction -= dx; // same as fraction -= 2*dx
                        }
                        x1 += stepx;
                        fraction += dy; // same as fraction -= 2*dy
                        putPixel(x1, y1);
                }
        } else {      // slope is not smaller than 45 degrees
                int fraction = dx - (dy >> 1);
                while (y1 != y2) {
                        if (fraction >= 0) {
                                x1 += stepx;
                                fraction -= dy;
                        }
                        y1 += stepy;
                        fraction += dx;
                        putPixel(x1, y1);
                }
        }
}

void Grafix::drawBox(int begin_x, int end_x, int begin_y, int end_y)
{
    lineBresenham(begin_x,begin_y,end_x,begin_y);
    lineBresenham(end_x,begin_y,end_x,end_y);
    lineBresenham(end_x,end_y,begin_x,end_y);
    lineBresenham(begin_x,end_y,begin_x,begin_y);
}

void Grafix::floatLine(float x1, float y1, float x2, float y2)
{
        // offset and stepping
        float dy = y2-y1;
        float dx = x2-x1;
        int stepx, stepy;

        int xe = (int)x2, ye = (int)y2;
        int x = (int)x1;
        int y = (int)y1;

        if (dy < 0) { stepy = -1; } else { stepy = 1; }
        if (dx < 0) { stepx = -1; } else { stepx = 1; }

        if (fabs(dx) > fabs(dy)) { // slope is less than 45 degrees
                float m = dy/dx;
                float b = y1-m*x1;

                while (x != xe) {
                        y = (int)(m*x+b);
                        putPixel(x, y);
                        x += stepx;
                }
                putPixel(xe, ye);
        } else {      // slope is not smaller than 45 degrees
                float m = dx/dy;
                float b = x1-m*y1;

                while (y != ye) {
                        x = (int)(m*y+b);
                        putPixel(x, y);
                        y += stepy;
                }
                putPixel(xe, ye);
        }
}

void Grafix::plotCircle(int x, int y, int x1, int y1)
{
        putPixel(x+x1, y+y1);
        putPixel(x-x1, y+y1);
        putPixel(x+x1, y-y1);
        putPixel(x-x1, y-y1);
        putPixel(x+y1, y+x1);
        putPixel(x-y1, y+x1);
        putPixel(x+y1, y-x1);
        putPixel(x-y1, y-x1);
}

void Grafix::circleBresenham(int x, int y, int radius)
{
        int x1 = 0;
        int y1 = radius;
        int p = 3 - 2*radius;

        while (x1 < y1) {
                plotCircle(x, y, x1, y1);
                if (p < 0)
                        p = p + 4*x1 + 6;
                else {
                        p = p + 4*(x1 - y1) + 10;
                        y1 = y1 - 1;
                }
                x1 = x1 + 1;
        }

        if ( x1 == y1 )
                plotCircle(x, y, x1, y1);
}

// ========================================================================
// PRIMITIVE DRAWING FUNCTIONS (WITH VIEWPORT TRANSFORMATION)
// ========================================================================

// draw a pixel on the canvas
void Grafix::drawPixel(double x, double y)
{
        x = vp_scale * x + vp_off_x;
        y = vp_scale * y + vp_off_y;

        putPixel((int)x, (int)y);
}

// draw a line on the canvas
void Grafix::drawLine(double x1, double y1, double x2, double y2)
{
        x1 = vp_scale * x1 + vp_off_x;
        y1 = vp_scale * y1 + vp_off_y;

        x2 = vp_scale * x2 + vp_off_x;
        y2 = vp_scale * y2 + vp_off_y;

        // use the bresenham algorithm for drawing a line
        lineBresenham((int)x1, (int)y1, (int)x2, (int)y2);
}

// draw a circle on the canvas
void Grafix::drawCircle(double x, double y, double radius)
{
        x = vp_scale * x + vp_off_x;
        y = vp_scale * y + vp_off_y;
        radius = radius * vp_scale;

        // draw actual circle
        circleBresenham((int)x, (int)y, (int)radius);
}

void Grafix::seedFillRecursive(int x, int y)
{
        int off = 0;
        int right, left;

        if ((x < 0) || (x >= resX) || (y < 0) || (y >= resY) ||
                    (canvas[y*resX+x] != ffBg)) return;

        while ((x+off >= 0) && (y >= 0) && (x+off < resX) && (y < resY) &&
               (canvas[y*resX+x+off] == ffBg)) putPixel(x+off++, y);

        right = off;
        off = -1;
        while ((x+off >= 0) && (y >= 0) && (x+off < resX) && (y < resY) &&
               (canvas[y*resX+x+off] == ffBg)) putPixel(x+off--, y);

        left = off;

        for(off = left+1; off < right; off++) {
                seedFillRecursive(x+off,y-1);
                seedFillRecursive(x+off,y+1);
        }
}

void Grafix::seedFill(double x, double y)
{
        int k = (int)(vp_scale * x + vp_off_x);
        int l = (int)(vp_scale * y + vp_off_y);

        ffBg = canvas[l*resX+k];
        if (ffBg == drawCol) return;
        seedFillRecursive(k, l);
}

// ========================================================================
// AUXILIARY FUNCTIONS
// ========================================================================

void Grafix::viewportTF(double x, double y, int &x_o, int &y_o)
{
        x_o = (int)(vp_scale * x + vp_off_x);
        y_o = (int)(vp_scale * y + vp_off_y);
}

void Grafix::inverseViewportTF(int x, int y, double &x_o, double &y_o)
{
        x_o = (x + 0.5 - vp_off_x) / vp_scale;
        y_o = (y + 0.5 - vp_off_y) / vp_scale;
}
