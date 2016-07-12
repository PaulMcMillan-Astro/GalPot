/*******************************************************************************
*                                                                              *
*  Potential.h                                                                 *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2006/07                                   *
* Lund Observatory, Lund University.                                           *
* address: Box 43, SE-221 00 Lund, Sweden                                      *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#ifndef _Pot_def_
#define _Pot_def_ 1

#include "Vector.h"
typedef Vector<double,3> Frequencies;


/*
   class Potential
   gives the base for user defined potentials; i.e. the user defines a class
   derived from Potential, which provides the operator() with 2 and 4 arguments
   (the latter for derivatives). An example is the logarithmic potential as 
   implemented in class LogPot.

   class DerPot
   is another base class derived from class Potential. Additionally to the 
   latter it provides the second derivatives of the potential.
*/

class Potential {
protected:
    double Lzsq;

public:

    void set_Lz    (const double Lz)    { Lzsq = Lz*Lz; }

    Potential      (const double Lz=0.) { set_Lz(Lz); }

    double Lzsquare() const             { return Lzsq; }

    virtual double operator()(				// returns Phi(R,z)
			      const double,		// given R
			      const double)const=0;	// and z

    virtual double operator()(				// returns Phi(R,z)
			      const double,		// given R
			      const double,		// and z
			      double&,			// also computes dPhi/dR
			      double&)const=0;		// and dPhi/dz

    virtual double RfromLc   (				// returns Rc,
			      const double,		// given Lz. possibly
			      double* =0)const=0;	// returns dRc/dLz.
    
    virtual double LfromRc   (				// returns Lc,
				const double,		// given R. possibly
				double* =0)const=0;	// returns dLz/dRc.
    
    virtual Frequencies KapNuOm(			// returns kappa,nu,Om
				const double)const=0;	// given R at z=0

    double eff(const double R, const double z) const
	{ if(Lzsq==0.) return (*this)(R,z);
          if(R==0.) {
	      cerr << " error in class Potential::eff: R=0 at non-zero Lz\n";
              exit(1); }
          return (*this)(R,z) + 0.5 * Lzsq/(R*R); }

    double eff(const double R, const double z, double& dPdR, double& dPdz) const
	{ 
	  if(Lzsq==0.) return (*this)(R,z,dPdR,dPdz);
    	  if(R==0.) {
              cerr << " error in class Potential::eff: R=0 at non-zero Lz\n";
              exit(1); }
    	  register double potential     = (*this)(R,z,dPdR,dPdz);
   	  register double Lzsq_over_Rsq = Lzsq/(R*R);
   	  dPdR                         -= Lzsq_over_Rsq / R;
    	  return potential + 0.5 * Lzsq_over_Rsq; }
};

#endif
