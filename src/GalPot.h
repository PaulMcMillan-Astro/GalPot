// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot.h                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// Copyright Walter Dehnen, 1996-2004                                          |
// e-mail:   walter.dehnen@astro.le.ac.uk                                      |
// address:  Department of Physics and Astronomy, University of Leicester      |
//           University Road, Leicester LE1 7RH, United Kingdom                |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Header file for class GalaxyPotential                                       |
//                                                                             |
// TO BE INCLUDED BY THE ENDUSER                                               |
//                                                                             |
// Version 0.0    15. July      1997                                           |
// Version 0.1    24. March     1998                                           |
// Version 0.2    22. September 1998                                           |
// Version 0.3    07. June      2001                                           |
// Version 0.4    22. April     2002                                           |
// Version 0.5    05. December  2002                                           |
// Version 0.6    05. February  2003                                           |
// Version 0.7    23. September 2004  fixed "find(): x out of range" error     |
// Version 0.8    04. February  2005  debugged error in GaussLegendre()        |
//                                    used WDutils::tupel                      |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// About the unit system                                                       |
// The unit system used throughout the classes and functions defined here and  |
// in included files is based on the following basic units                     |
// 	unit of length:		1 kilo parsec (kpc)                            |
//	unit of time:		1 mega year   (Myr)                            |
// 	unit of mass:		1 solar mass  (Msun)                           |
// This implies the following dimensions                                       |
//                                                                             |
// quantity        dimension / seize                using other units          |
//-----------------------------------------------------------------------------+
// angular vel.  1 Myr^-1	                   = 977.775320024919 km/s/kpc |
// velocity      1 kpc/Myr                         = 977.775320024919 km/s     |
// action/mass   1 kpc^2/Myr                       = 977.775320024919 kpc*km/s |
// potential     1 kpc^2/Myr^2                     = 956044.576449833 (km/s)^2 |
// acceleration  1 kpc/Myr^2                                                   |
// G             4.49865897 E-12 kpc^3/Myr^2/Msun                              |
// 4 Pi G        5.65318158 E-11 kpc^3/Myr^2/Msun                              |
//-----------------------------------------------------------------------------+
#ifndef GalPot_h
#define GalPot_h

#include <iomanip>
#include <iostream>
#include <cmath>                                 // v0.4                        
#include <cstdlib>                               // v0.4
#include <string>
#include "Units.h"
#include "Vector.h"
#include "Potential.h"

  //----------------------------------------------------------------------------
  // define the maximum l used in the multipole expansion                       
  // After adjustment re-compile GalPot.cc for the change to take effect        
  //----------------------------------------------------------------------------

#ifdef GalPot_cc
const int LMAX=80;		// maximum l for the multipole expansion        
#endif

//----------------------------------------------------------------------------
// Set some values used as default in the constructors below                  
//----------------------------------------------------------------------------

const int    NRAD=201;	      // DEFAULT number of radial points in Multipole 
const double RMIN=1.e-4*Units::kpc,// DEFAULT min radius of logarithmic radial grid
             RMAX=1.e3*Units::kpc; // DEFAULT max radius of logarithmic radial grid
//----------------------------------------------------------------------------
// Include the definitions for all auxiliary functions. These contain a long  
// tail of dependencies needed. It seems unavoidable to have them all visible 
// to the enduser, if GalaxyPotential below is defined as a proper class, e.g.
// the construction of not just one object is possible etc .                  
//                                                                            
// Most important for the enduser is the meaning of DiskPar and SphrPar.      
// These are Vectors of 5 and 6 doubles, respectively, holding the parameters 
// for one disk or spheroid component. The meaning of them is as follows      
// DiskPar[0]   is the surface density normalisation Sigma_0 [Msun/kpc^2]     
// DiskPar[1]   is the scale length R_d [kpc]                                 
// DiskPar[2]   is the scale height h [kpc]. For h<0 an isothermal (sech^2)   
//		  profile is used, for h>0 an exponential one, and for h=0 the  
//		  disk is infinitesimal thin.                                   
// DiskPar[3]   is the inner cut-off radius R_m [kpc]                         
// DiskPar[4]   is eps. A term eps*cos(pi*R/R_d) is added to the exponent.       
//                                                                            
// SphrPar[0]   is the density normalization rho_0 [Msun/kpc^3]               
// SphrPar[1]   is the axis ration q                                          
// SphrPar[2]   is the inner power slope gamma                                
// SphrPar[3]   is the outer power slope beta                                 
// SphrPar[4]   is the transition radius r_0 [kpc]                            
// SphrPar[5]   is the outer cut-off radius r_t [kpc]                         
//----------------------------------------------------------------------------

typedef Vector<double,5> DiskPar;  
typedef Vector<double,6> SphrPar;
class PotResidual {
public:
  // density at given (R,z)
  virtual double Density  (const double, const double)                const=0;
  // residual density (input for multipole expansion) at given (r,sin/cos(th))
  virtual double Residual (const double, const double, const double)  const=0;
};
//////////////////////////////////////////////////////////////////////////////
class DiskAnsatz : public PotResidual {
private:
  double               S0, Rd, zd, R0, eps;           // defining  variables  
  int                  thin, hollow, isothermal;      // auxiliary variable   
  double               Rd2, zdoRd, fac, R0oRd;        // auxiliary variables  
  double               mass_integrand(const double) const;
public:
  void                 setup(const DiskPar&);
  DiskAnsatz           () {}
  DiskAnsatz           (const DiskPar&d) { setup(d); }
  // return potential (and its gradient) of part not in residual
  double operator()    (const double, const double, const double, double* =0)
    const;
  double      Laplace       (const double, const double) const;
  Frequencies      kapnuom       (const double) const;
  bool        is_thin       () const { return thin; }
  bool        is_hollow     () const { return hollow; }
  double      mass          (const double=0.) const;
  double      SurfaceDensity(const double) const;
  double      Density       (const double, const double) const;
  double      Residual      (const double, const double, const double) const;
  void        DescribePot   (ostream&) const;
  DiskPar     parameter     () const;
};
inline DiskPar DiskAnsatz::parameter() const
{
  DiskPar p;
  p[0]=S0; p[1]=Rd; p[2]=(isothermal)? -zd : zd; p[3]=R0; p[4]=eps;
  return p;
}
inline void DiskAnsatz::DescribePot(ostream& to) const
{
    to<< "potential due to mass density:\n "
      <<std::setprecision(3)<<S0<<Units::mass_unit<<'/'<<Units::length_unit<<"^2 Exp{-";
    if(hollow) 
	to<<std::setprecision(3)<<R0<<Units::length_unit<<"/R-";
    to<<"R/"<<std::setprecision(3)<<Rd<<Units::length_unit;
    if(isothermal)
	to<<"} sech^2{z/"
	  <<std::setprecision(3)<<(2*zd)<<Units::length_unit<<"}/"
	  <<std::setprecision(3)<<(4*zd)<<Units::length_unit;
    else if(!thin)
	to<<"-|z|/"<<std::setprecision(3)<<zd<<Units::length_unit<<"}/"
	  <<std::setprecision(3)<<(2*zd)<<Units::length_unit;
    to<<'\n';
}
//////////////////////////////////////////////////////////////////////////////
class Disks {
protected:
  int                   nd;
  DiskAnsatz            *D, *Dup;
  void reset            (const int, const DiskPar*);
public:
  Disks                 (std::istream&);
  Disks                 (const Disks&); 
  Disks                 (const int, const DiskPar*);
  ~Disks                 () { delete[] D; }
  bool    all_hollow    () const;
  bool    none_hollow   () const;
  double  Mass          (const double=0.) const;
  double  SurfaceDensity(const double) const;
  double  operator()    (const double, const double) const;
  double  operator()    (const double, const double, double&, double&) const;
  double  Laplace       (const double, const double) const;
  double  Density       (const double, const double) const;
  double  Residual      (const double, const double, const double) const;
  int     NumberofDisks () const { return nd; }
  DiskPar Parameter     (const int i) const { return (D+i)->parameter(); }
  void DescribePot(ostream& ) const;
};
inline bool Disks::none_hollow() const { 
  if(nd==0) return 1;
  for(register DiskAnsatz *p=D; p<Dup; p++) if( (p->is_hollow()) ) return 0;
  return 1;
}
inline bool Disks::all_hollow() const {
  if(nd==0) return 0;
  for(register DiskAnsatz *p=D; p<Dup; p++) if( !(p->is_hollow()) ) return 0;
  return 1;
}
inline double Disks::Mass(const double r) const { 
  if(nd==0) return 0.; 
  register double R=0.;
  for(register DiskAnsatz *p=D; p<Dup; p++) R += p->mass(r);
  return R;
}
inline double Disks::SurfaceDensity(const double a) const {
  if(nd==0) return 0.;
  register double R=0.;
  for(register DiskAnsatz *p=D; p<Dup; p++) R += p->SurfaceDensity(a);
  return R;
}
inline double Disks::operator() (const double R, const double z) const {
  if(nd==0) return 0.;
  register double r=hypot(R,z), pot=0.;
  for(register DiskAnsatz *p=D; p<Dup; p++) pot+= (*p)(R,z,r);
  return pot;
}
inline double Disks::operator() (const double R, const double z,
				 double& dR, double& dz) const {
  if(nd==0) {
    dR = 0.;
    dz = 0.;
    return 0.;
  }
  register double d[2], r=hypot(R,z), pot=0.;
  register DiskAnsatz *p=D;
  for(dR=dz=0.; p<Dup; p++)
    { pot += (*p)(R,z,r,d); dR+=d[0]; dz+=d[1]; }
  return pot;
}
inline double Disks::Laplace(const double a, const double b) const { 
  if(nd==0) return 0.;
  register double L=0.;
  for(register DiskAnsatz *p=D; p<Dup; p++) L += p->Laplace(a,b);
  return L;
}
inline double Disks::Density(const double a, const double b) const {
  if(nd==0) return 0.;
  register double R=0.;
  for(register DiskAnsatz *p=D; p<Dup; p++) R += p->Density(a,b);
  return R;
}
inline double Disks::Residual(const double a, const double b, const double c)
  const {
  if(nd==0) return 0.;
  register double R=0.;
  for(register DiskAnsatz *p=D; p<Dup; p++) R += p->Residual(a,b,c);
  return R;
}
inline void   Disks::DescribePot(ostream& to) const
{
    register DiskAnsatz *p=D;
    for(; p<Dup; p++) p->DescribePot(to);
}
//////////////////////////////////////////////////////////////////////////////
class SpheroidDensity : public PotResidual {
private:
  double rh0, q, gam, bet, r0, rcut;              // defining variables       
  double beg, qi, r0i, rci;                       // auxiliary variables      
  double mass_integrand(const double) const;      // auxiliary function       
public:
  void    setup(const SphrPar&);
  SpheroidDensity        () {}
  SpheroidDensity        (const SphrPar &s)       { setup(s); }
  bool    cut_off        () const { return (rci>0.);  }
  double  scale_density  () const { return rh0;  }
  double  inner_power    () const { return gam;  }
  double  outer_power    () const { return (rci)? 1.e3 : bet; }
  double  mass           (const double) const;
  double  Density        (const double, const double) const;
  double  Residual       (const double, const double, const double) const;
  SphrPar parameter      () const;
  void   DescribePot     (ostream&) const;
};
inline double SpheroidDensity::Residual(const double r, const double st,
					const double ct) const {
  return Units::fPiG * Density(r*st,r*ct);
}
inline SphrPar SpheroidDensity::parameter() const {
  SphrPar p;
  p[0]=rh0;
  p[1]=q;
  p[2]=gam;
  p[3]=bet;
  p[4]=r0;
  p[5]=rcut;
  return p;
}
inline void SpheroidDensity::DescribePot(ostream& to) const
{
    to<<" Spheroid: rh0="<<rh0
	       <<", c/a="<<q
	       <<", gamma="<<gam
	       <<", beta="<<bet
	       <<", r0="<<r0;
    if(rci) to<<", rcut="<<rcut;
    to<<'\n';
}

//////////////////////////////////////////////////////////////////////////////
class Spheroids {
  int                      ns;
  SpheroidDensity          *S, *Sup;
protected:
  void      reset          (const int, const SphrPar*);
public:
  Spheroids                (std::istream&);
  Spheroids                (const int, const SphrPar*);
  Spheroids                (const Spheroids&);
  Spheroids& operator=     (const Spheroids&);
  ~Spheroids                () { delete[] S; }
  bool    massive          () const;
  double  beta             () const;
  double  gamma            () const;
  double  Mass             (const double) const;
  double  Density          (const double, const double) const;
  double  Residual         (const double, const double, const double)const;
  int     NumberofSpheroids() const { return ns; }
  SphrPar Parameter        (const int i) const { return (S+i)->parameter(); }
  void DescribePot(ostream& ) const;
};
inline bool Spheroids::massive() const {
  if(ns==0) return 0;
  for(register SpheroidDensity *p=S; p<Sup;p++)
    if(p->scale_density()) return true;
  return false;
}
inline double Spheroids::Mass(const double a) const {
  if(ns==0.) return 0.;
  register double R=0.;
  for(register SpheroidDensity *p=S; p<Sup; p++) R += p->mass(a);
  return R;
}
inline double Spheroids::Density(const double a, const double b) const {
  if(ns==0.) return 0.;
  register double R=0.;
  for(register SpheroidDensity *p=S; p<Sup; p++) R += p->Density(a,b);
  return R;
}
inline double Spheroids::Residual(const double a, const double b,
				  const double c) const {
  if(ns==0) return 0.;
  register double R=0.;
  for(register SpheroidDensity *p=S; p<Sup; p++) R += p->Residual(a,b,c);
  return R;
}
inline void   Spheroids::DescribePot(ostream& to) const
{
    register SpheroidDensity *p=S;
    for(; p<Sup; p++) p->DescribePot(to);
}

//////////////////////////////////////////////////////////////////////////////
class Multipole {
private:
  int         LR,K[2];
  double      Rmin, Rmax, gamma, beta, Phi0;
  double      lRmin, lRmax, g2;
  double      lzmin, lzmax, tg3, g3h;
  double      *logr, *lLc, *d2R, *d2L; 
  double      *X[2], **Y[3], **Z[4];
  void        AllocArrays();
  void        setup(const double, const double,       // r_min, r_max         
		    const double, const double,       // gamma, beta          
		    PotResidual const*);              // providing rho(x)     
public:
  Multipole (const int,                               // points on log grid   
	     const double, const double,              // r_min, r_max         
	     const double, const double,              // gamma, beta          
	     PotResidual const*,                      // providing rho(x)     
	     const int =1);                           // routines LfromRc...  
  void reset(const double, const double,              // r_min, r_max         
	     const double, const double,              // gamma, beta          
	     PotResidual const*,                      // providing rho(x)     
	     const int =1);                           // routines LfromRc...  
  ~Multipole();
  double      operator()(const double,const double,const double, double* =0) 
    const;
  double      vcsquare  (const double)		const;
  double      vcsquare  (const double, double&)	const;
  double      RfromLc   (const double, double* =0)	const;
  double      LfromRc   (const double, double* =0) 	const;
  double      Laplace   (const double, const double)  const;
  Frequencies      kapnuom   (const double)		const;
};
inline double Multipole::vcsquare(const double R) const
{ 
  double d[2];
  operator() (R,0.,1.,d);
  return R*d[0];
}




class GalaxyPotential : public  PotResidual,
			public  Potential,
			private Disks,
			private Spheroids {
private:
  
  Multipole M;
  
  double Residual      (const double, const double, const double) const;
  // this function is not intended for the enduser, it is needed in the       
  // construction of GalaxyPotential itself                                   
  
public:
  
  virtual ~GalaxyPotential() {}
  
  //--------------------------------------------------------------------------
  // constructors and related functions                                       
  //--------------------------------------------------------------------------
  GalaxyPotential(std::istream&);
  // constructor from input stream. Use this constructor to establish the     
  // potential of one of the models as shown in the following code fragment.  
  //                                                                          
  // ifstream from("Model.pot"); 	// file `Model.pot' contains the data   
  // GalaxyPotenial Phi(from);	// read from file and construct object  
  // from.close();			// close file                           

  //  GalaxyPotential(std::string&);
  // constructor from file with name given in string.

  
  GalaxyPotential(const int, const DiskPar*,	// No & parameters of disks     
		  const int, const SphrPar*,	// No & parameters of spheroids 
		  const double = RMIN,	// min radius of radial grid    
		  const double = RMAX, 	// max radius of radial grid    
		  const int    = NRAD);	// No of radial grid points     
  // constructor from parameters. See above for the meaning of the parameters 
  
  void   reset   (const int, const DiskPar*,	// No & parameters of disks     
		  const int, const SphrPar*, 	// No & parameters of spheroids 
		  const double = RMIN,	// min radius of radial grid    
		  const double = RMAX);	// max radius of radial grid    
  // resets to new parameters (arguments as the above constructor), the       
  // number of radial grid points remains fixed at the old value              
  
  //--------------------------------------------------------------------------
  // applications (are all const member functions)                            
  //--------------------------------------------------------------------------
  // information on the disks alone                                           
  //--------------------------------------------------------------------------
  
  double DisksDensity(const double, const double) const;
  // returns the disks' volume density at some (R,z) (1st & 2nd argument)     
  
  double DisksSurfaceDensity(const double) const;
  // returns the disks' surface density at some radius (1st argument)         
  // = density integrated over z from -oo to oo.                              
  
  double DisksMass(const double=0.) const;
  // returns the disks' mass inside some radius (1st argument)                
  // = 2*Pi*R*density integrated over z from -oo to oo, and R from 0 to R.    
  // if called with no or zero argument, the total mass is returned           
  
  int NumberofDisks() const;
  // returns the number of disk components                                    
  
  DiskPar DiskParameter(const int) const;
  // returns the parameters of the ith (1st argument) disk                    
  
  //--------------------------------------------------------------------------
  // information on the spheroids alone                                       
  //--------------------------------------------------------------------------
  
  double SpheroidsDensity(const double, const double) const;
  // returns the spheroids' volume density at some (R,z) (1st & 2nd argument) 
  
  double SpheroidsMass(const double) const;
  // returns the spheroids' total mass inside some radius (1st argument)      
  // = 4*Pi*q*m^2*density integrated over m from 0 to R                       
  
  int NumberofSpheroids() const;
  // returns the number of spheroid components                                
  
  SphrPar SpheroidParameter(const int) const;
  // returns the parameters of the ith (1st argument) spheroid                
  
  //--------------------------------------------------------------------------
  // information on the total                                                 
  //--------------------------------------------------------------------------
  
  double Density(const double, const double) const;
  // returns the (input) density in Msun/kpc^3 given (R,z) in kpc.            
  
  double vcsquare(const double) const;
  // return the circular speed squared in (kpc/Myr)^2 at R given in kpc and   
  // z=0.  vc^2 is defined by R*dPhi/dR, which can become negative, e.g. in   
  // the central parts of a hollow disk                                       
  
  double Mass(const double) const;
  // returns the total mass inside some radius (1st argument)                 
  // this is simply the sum of DisksMass() and SpheroidsMass() above          

  void OortConstants(const double, double&, double&) const;
  // given R in kpc (1st argument) returns Oort's constants A, B as 2nd and   
  // 3rd argument. Units of the latter are 1/Myr                              
  
  double operator() (const double, const double) const;
  // returns Phi in (kpc/Myr)^2 at (R,z) given in kpc                         

  double operator() (const double, const double, double&, double&) const;
  // returns Phi in (kpc/Myr)^2 at (R,z) given in kpc.                        
  // additionally, on return the 3rd and 4th argument contain the derivatives 
  // dPhi/dR and dPhi/dz in units of kpc/Myr^2                                
  
  double RfromLc       (const double, double* =0) const;
  double LfromRc       (const double, double* =0) const;
  
  double Laplace(const double, const double) const;
  // returns Laplace(Phi) in 1/Myr^2 given (R,z) in kpc. This is not          
  // necessarily identical to 4 Pi G times GalaxyPot::Density() above, as the 
  // latter returns the input density, while this routine gives the Laplace   
  // of the potential as evaluated. However, both numbers should agree in the 
  // limit of infinite many radial grid points, multipoles, and infinite      
  // numerical accuracy                                                       
  
  Frequencies KapNuOm  (const double) const;
  // given R, the epicycle frequencies kappa (radial), nu (vertical), and     
  // omega (azimuthal) for the circular orbit through (R,z=0) are returned.   
  // units for the frequencies are 1/Myr                                      
  void DescribePot(ostream& ) const;
};

//////////////////////////////////////////////////////////////////////////////
// Constructors and related (all inline)                                      
//////////////////////////////////////////////////////////////////////////////
inline double GalaxyPotential::Residual(const double a,
					const double b,
					const double c) const {
  return Disks::Residual(a,b,c) + Spheroids::Residual(a,b,c);
}
//----------------------------------------------------------------------------
inline GalaxyPotential::GalaxyPotential(std::istream &from)
  : Disks(from),
    Spheroids(from),
    M(NRAD,RMIN,RMAX,Spheroids::gamma(),Spheroids::beta(),this)
{}
//----------------------------------------------------------------------------
inline GalaxyPotential::GalaxyPotential(const int Nd, const DiskPar* pd,
					const int Ns, const SphrPar* ps,
					const double rmin, const double rmax,
					const int k)
  : Disks(Nd,pd),
    Spheroids(Ns,ps),
    M(k,rmin,rmax,Spheroids::gamma(),Spheroids::beta(),this)
{}
//----------------------------------------------------------------------------
inline void GalaxyPotential::reset(const int Nd, const DiskPar* pd,
				   const int Ns, const SphrPar* ps,
				   const double rmin, const double rmax) {
  Disks::reset(Nd,pd);
  Spheroids::reset(Ns,ps);
  M.reset(rmin,rmax,Spheroids::gamma(),Spheroids::beta(),this);
}
//////////////////////////////////////////////////////////////////////////////
// information on the disks alone (all inline)                                
//////////////////////////////////////////////////////////////////////////////
inline double GalaxyPotential::DisksDensity(const double R, const double z)
  const {
  return Disks::Density(R,z);
}
//----------------------------------------------------------------------------
  inline double GalaxyPotential::DisksSurfaceDensity(const double R) const {
    return Disks::SurfaceDensity(R);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::DisksMass(const double R) const {
    return Disks::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline int GalaxyPotential::NumberofDisks() const {
    return Disks::NumberofDisks();
  }
  //----------------------------------------------------------------------------
  inline DiskPar GalaxyPotential::DiskParameter(const int i) const {
    return Disks::Parameter(i);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the spheroids alone (all inline)                            
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::SpheroidsDensity(const double R, const double z)
    const {
    return Spheroids::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::SpheroidsMass(const double R) const {
    return Spheroids::Mass(R);
  }
  //----------------------------------------------------------------------------
  inline int GalaxyPotential::NumberofSpheroids() const {
    return Spheroids::NumberofSpheroids();
  }
  //----------------------------------------------------------------------------
  inline SphrPar GalaxyPotential::SpheroidParameter(const int i) const {
    return Spheroids::Parameter(i);
  }
  //////////////////////////////////////////////////////////////////////////////
  // information on the total (most are non-inline functions)                   
  //////////////////////////////////////////////////////////////////////////////
  inline double GalaxyPotential::Density(const double R, const double z) const {
    return Disks::Density(R,z) + Spheroids::Density(R,z);
  }
  //----------------------------------------------------------------------------
  inline double GalaxyPotential::Mass(const double R) const {
    return Disks::Mass(R) + Spheroids::Mass(R);
  }
  //----------------------------------------------------------------------------
inline double GalaxyPotential::vcsquare(const double R) const {
    return M.vcsquare(R);
}

inline double GalaxyPotential::RfromLc(const double L, double* dR) const
{
  return M.RfromLc(fabs(L),dR);
}

inline double GalaxyPotential::LfromRc(const double R, double* dR) const
{
    return M.LfromRc(R,dR);
}
inline void GalaxyPotential::DescribePot(ostream& to) const
{
    Disks::DescribePot(to);
    Spheroids::DescribePot(to);
}


#endif  // #ifndef GalPot_h
