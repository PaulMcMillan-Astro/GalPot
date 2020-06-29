/***************************************************************************//**
\file MiyamotoNagaiPot.h
\brief Contains class MiyamotoNagaiPotential
Which is a Miyamoto-Nagai potential. The clue is in the name.

*                                                                              *
*  MiyamotoNagaiPot.h                                                          *
*                                                                              *
* C++ code written by Paul McMillan, 2006-                                  
*                                                                              *
*******************************************************************************/

#ifndef _MiyNag_h_
#define _MiyNag_h_ 1

#include <iostream>
#include "Potential.h"
#include "Units.h"

/** \brief A Miyamoto-Nagai potential.

Input values are M, a, b. See Binney & Tremaine (2008), eq 2.69 for
the potential.

 */
class MiyamotoNagaiPotential : public Potential {
public:
  const double A,ABoB,Bq,ABq,GM;

  //public:
    MiyamotoNagaiPotential(const double=1., const double=1.,
	       const double=1.);
    double operator() (const double, const double) const;
    double operator() (const double, const double, double&, double&) const;
    //double RfromLc(const double, double* = 0) const {return 0.;} // Not implimented yet
    //double LfromRc(const double, double* = 0) const {return 0.;} // Not implimented yet
    Frequencies KapNuOm(const double) const;
    friend ostream& operator<< (ostream&, const MiyamotoNagaiPotential&);
};

inline MiyamotoNagaiPotential::MiyamotoNagaiPotential(const double m, const double a,
				    const double b)
		     : GM(Units::G*m), A(a), ABoB((A+b)/b), Bq(b*b),
		     ABq((A+b)*(A+b))
{ }


#endif
