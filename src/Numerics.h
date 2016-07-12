/*******************************************************************************
*                                                                              *
* Numerics.h                                                                   *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
* Lund Observatory, Lund University.                                           *
* address: 1Keble Road, Oxford, OX1 3NP, United Kingdom.                       *
* e-mail:  dehnen@thphys.ox.ac.uk                                              *
*                                                                              *
*******************************************************************************/

#ifndef _Numerics_def_
#define _Numerics_def_ 1

#include "Vector.h"
#include "Matrix.h"

#include "Pi.h"
#include "WDMath.h"
#include "Inline.h"
#include <algorithm>
using std::max;
using std::min;
#include "Numerics.templates"

////////////////////////////////////////////////////////////////////////////////
// here only the non-inline non-template functions are listed: /////////////////
////////////////////////////////////////////////////////////////////////////////


void GaussLegendre(double*,double*,const int);
void LegendrePeven(double*,const double,const int);
void dLegendrePeven(double*,double*,const double,const int);

double qbulir(double(*)(double),const double,const double,const double,double&);
inline double qbulir(double(*func)(double),const double a,const double b,
		     const double eps)
{ 
    double err; 
    return qbulir(func,a,b,eps,err);
}

#endif
