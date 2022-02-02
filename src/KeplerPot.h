/***************************************************************************//**
\file KeplerPot.h
\brief Contains class KeplerPotential.

Gives the Kepler (point mass) potential.

*                                                                              *
* KeplerPot.h                                                                   *
*                                                                              *
* C++ code written by Paul McMillan, 2014-
*                                                                              *
*******************************************************************************/

#ifndef KeplerPotential_GalPot_
#define KeplerPotential_GalPot_ 1

#include <iostream>
#include <cmath>
#include "Potential.h"
#include "Units.h"

/** \brief A Kepler potential.

    Input parameter is mass.

 */

class KeplerPotential : public Potential {
  double GM;
public:
  KeplerPotential();
  ~KeplerPotential() {;}
  KeplerPotential(double);
  double operator() (const double, const double) const;
  double operator() (const double, const double, double&, double&) const;
  //double RfromLc(const double, double* = 0) const {return 0.;} // Not implimented yet
  //double LfromRc(const double, double* = 0) const  {return 0.;} // Not implimented yet
  Frequencies KapNuOm(const double) const;
};


inline KeplerPotential::KeplerPotential(double M) :
		      GM(M*Units::G)

{

}

#endif
