//
// C++ Implementation: hpmath
//
// Description: 
//
//
// Author: Jan Fostier <jan.fostier@intec.ugent.be>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "hpmath.h"
#include <cmath>
#include <cfloat>
#include <cassert>

double omxn(double x, int n)
{
    // a) we can evaluate 1-x without loss of precision
    if (std::abs(x) > 0.1) return pow(1.0-x,n);

    // b) we cannot => use taylor expansion
    const double epsilon = std::abs(x*DBL_EPSILON);

    double term=1.0;
    double result=term;

    int i=0;
    while ((std::abs(term) > epsilon)  and i!=n) {
        i++;
        term=-term*(n-i+1.0)/(i*1.0)*x;
        result+=term;
    }
    return result;
}

double logopx(double x)
{
	// a) we can evaluate x+1 without loss of precision
	if (std::abs(x) > 0.1) return log(1.0+x);

	// b) we cannot evaluate x+1 without loss of precision
	// We use log(1+x) = x - x^2/2 + x^3/3 - x^4/4 + ...
	const double epsilon = std::abs(x*DBL_EPSILON);
	double xpn = x;
	double term = xpn;
	double result = term;
	int n = 2;

	while (std::abs(term) > epsilon) {
		xpn *= -x;
		term = xpn/(n++);
		result += term;
	}

	return result;
}

double ompowopxn(double x, int n)
{
	// this function does not support negative n
	assert(n >= 0);

	// a) if we can evaluate x+1 without loss of precision
	if (std::abs(x) > 0.1) return 1.0-pow(1.0+x, n);

	// b) if we cannot calculate the binomial coefficients
	if (n > 30) return 1.0-pow(1.0+x, n);

	// b) we cannot evaluate x+1 without loss of precision
	int binomialCoef = n;
	double xpn = x;
	double term = binomialCoef*xpn;	
	double result = term;

	for (int i = n-2, j = 2; j <= n; j++, i--) {
		binomialCoef *= (i+1);
		binomialCoef /= j;
		xpn *= x;

		term = binomialCoef*xpn;
		result += term;
	}

	return -result;
}

double omexpx(double x)
{
	// a) if we can evaluate 1-exp(x) without loss of precision
	if (std::abs(x) > 0.1) return 1.0-exp(x);

	// b) we cannot evaluate 1-exp(x) without loss of precision
	// We use 1-exp(x) = -x - x^2/2! - x^3/3! - x^4/4! + ...
	const double epsilon = std::abs(x*DBL_EPSILON);
	double xpn = x;
	double term = xpn;
	double factor = 1.0;
	double result = term;
	int n = 2;

	while (std::abs(term) > epsilon) {
		xpn *= x;
		factor /= (n++);

		term = factor*xpn;
		result += term;
	}

	return -result;
}
