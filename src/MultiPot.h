/***************************************************************************//**
\file MultiPot.h
\brief Contains class MultiPotential.

Combines multiple potentials into one.

*                                                                              *
* MultiPot.h                                                                   *
*                                                                              *
* C++ code written by Paul McMillan, 2014-
*                                                                              *
*******************************************************************************/

#ifndef MultiPotential_GalPot_
#define MultiPotential_GalPot_ 1

#include <iostream>
#include <cmath>
#include "Potential.h"

/** \brief A sum of multiple potentials. The potentials can be defined separately.

    Input parameter is a pointer to an array of Potential pointers.

 */

class MultiPotential : public Potential {
  int npot;
  Potential **PotList;
public:
  MultiPotential();
  ~MultiPotential() {;}
  MultiPotential(Potential**, const int);
  double operator() (const double, const double) const;
  double operator() (const double, const double, double&, double&) const;
  //double RfromLc(const double, double* = 0) const;
  //double LfromRc(const double, double* = 0) const;
  Frequencies KapNuOm(const double) const;
};


inline MultiPotential::MultiPotential(Potential **inPotList,
				      const int np) :
		      PotList(inPotList), npot(np)

{

}

#endif
