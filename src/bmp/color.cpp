#include "color.h"
#include "stdlib.h"
#include "cstdlib"
#include <iostream>
using namespace std;
// ========================================================================
// COLOR
// ========================================================================

bool Color::operator==(const Color &rhs) const
{
        if (r != rhs.r) return false;
        if (b != rhs.b) return false;
        if (g != rhs.g) return false;

        return true;
}

bool Color::operator!=(const Color &rhs) const
{
        if (r != rhs.r) return true;
        if (b != rhs.b) return true;
        if (g != rhs.g) return true;

        return false;
}

void Color::getIntegerRGB(uint& red, uint& green, uint& blue) const
{
    red = (uint)r;
    green =(uint)g;
    blue = (uint)b;
}

void Color::getRandomColor()
{
    int rn1,rn2,rn3;
    bool badColor=false;
    do {
        rn1=rand()%256;
        rn2=rand()%256;
        rn3=rand()%256;

        r=(unsigned char) rn1;
        g=(unsigned char) rn2;
        b=(unsigned char) rn3;

        //if random color is close to white or black restart
        if ((rn1+rn2+rn3) > (3*255-100))
            badColor=true;
        else if ((rn1+rn2+rn3) < 100)
            badColor=true;
        else
            badColor=false;

    } while (badColor);
}
