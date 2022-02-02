/*******************************************************************************
*                                                                              *
*  Potential.h                                                                 *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2006/07                                   *
*                                                                              *
*******************************************************************************/

#ifndef Pot_GalPot_
#define Pot_GalPot_ 1

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
			      double* =0)const;	// returns dRc/dLz.

    virtual double LfromRc   (				// returns Lc,
  		      const double R,		// given R. possibly
  		      double* dR=0) const {        	// returns dLz/dRc.
      double dPR,dPz,P;
      P = (*this)(R,0.,dPR,dPz);
      return sqrt(R*R*R*dPR);
    }
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
    	  double potential     = (*this)(R,z,dPdR,dPdz);
   	  double Lzsq_over_Rsq = Lzsq/(R*R);
   	  dPdR                         -= Lzsq_over_Rsq / R;
    	  return potential + 0.5 * Lzsq_over_Rsq; }
};


inline double Potential::RfromLc(const double L, double* dR) const
{
  bool more=false;
  double R,lR=0.,dlR=0.001,z,dPR,dPz,P,LcR,oldL;
  R=exp(lR);
  P= (*this)(R,0.,dPR,dPz);
  LcR=pow(R*R*R*dPR,0.5);
  if(LcR == L) return R;
  if(L>LcR) more=true;
  oldL=LcR;
  for( ; ; ){
    lR += (more)? dlR : -dlR;
    R=exp(lR);
    P= (*this)(R,0.,dPR,dPz);
    LcR=pow(R*R*R*dPR,0.5);
    if(LcR == L) return R;
    if((L< LcR && L>oldL) ||(L>LcR && L<oldL)){
	R=(more)? exp(lR-0.5*dlR) : exp(lR+0.5*dlR);
	return R;}
    oldL=LcR;
  }

}


#endif
