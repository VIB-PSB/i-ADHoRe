#ifndef COLOR_H
#define COLOR_H

// ========================================================================
// COLOR
// ========================================================================

// color class
class Color{

public:
        Color() : r(0), g(0), b(0) {};
        Color(unsigned char r_, unsigned char g_, unsigned char b_) :
                        r(r_), g(g_), b(b_) {};

        /**
         * Check whether two colors equal
         * @param rhs Right hand side color
         * @return True of false
         */
        bool operator==(const Color &rhs) const;

        /**
         * Check whether two colors differ
         * @param rhs Right hand side color
         * @return True of false
         */
        bool operator!=(const Color &rhs) const;

        /**
         * Get the blue component
         * @return Blue component
         */
        unsigned char getBlue() { return b; }

        /**
         * Get the green component
         * @return Green component
         */
        unsigned char getGreen() { return g; }

        /**
         * Get the red component
         * @return Red component
         */
        unsigned char getRed() { return r; }

        void getIntegerRGB(unsigned int& red, unsigned int& green, unsigned int& blue) const;

        /**
        * Converts this to a random color (3 times rand%256)
        *NOTE randomcolor cannot be too close to black or white!
        */
        void getRandomColor();

private:
        unsigned char r, g, b;
};

// ========================================================================
// CONSTANTS
// ========================================================================

const Color red(0xff, 0 , 0);
const Color green(0, 0xff, 0);
const Color blue(0, 0, 0xff);
const Color yellow(0xff, 0xff, 0);
const Color cyan(0, 0xff, 0xff);
const Color purple(0xff, 0, 0xff);
const Color gray(0x90, 0x90, 0x90);
const Color white(0xff, 0xff, 0xff);
const Color black(0, 0, 0);

#endif
