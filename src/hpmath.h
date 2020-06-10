//
// C++ Interface: hpmath
//
// Description: 
//
//
// Author: Jan Fostier <jan.fostier@intec.ugent.be>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef HPMATH_H
#define HPMATH_H


double omxn(double x, int n);
/**
 * Evaluates (1-x)^n in a numerically safe manner, even for tiny x
 * @param x Argument x
 * @return (1-x)^n
 */
double omxn(double x, int n);

/**
 * Evaluates log(1+x) in a numerically safe manner, even for tiny x
 * @param x Argument x
 * @return Natural logarithm of 1+x
 */
double logopx(double x);

/**
 * Evaluates 1-(1+x)^n in a numerically safe manner, even for tiny x
 * @param x Argument x
 * @param n Exponential with n >= 0
 * @return Evaluation of 1-(1+x)^n
 */
double ompowopxn(double x, int n);

/**
 * Evaluates 1-exp(x) in a numerically safe manner, even for tiny x
 * @param x Argument x
 * @return Evaluation of 1-exp(x)
 */
double omexpx(double x);

#endif
