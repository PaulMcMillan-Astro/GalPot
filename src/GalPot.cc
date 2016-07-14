// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalPot.cc                                                                   |
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
// source code for class GalaxyPotential and dependencies                      |
//                                                                             |
// TO BE COMPILED AND LINKED TO THE ENDUSER                                    |
//                                                                             |
// Version 0.0    15. July      1997                                           |
// Version 0.1    24. March     1998                                           |
// Version 0.2    22. September 1998                                           |
// Version 0.3    07. June      2001                                           |
// Version 0.4    22. April     2002                                           |
// Version 0.5    05. December  2002                                           |
// Version 0.6    05. February  2003                                           |
// Version 0.7    23. September 2004  fixed "find(): x out of range" error     |
// Version 0.8    24. June      2005  explicit construction of tupel           |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef GalPot_cc
#define GalPot_cc
#include "GalPot.h" 
#include "FreeMemory.h"
//#include "Numerics.h"
#include "Pspline.h"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// class DiskAnsatz                                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//   rho(R,z) = Sig(R) * h(z)                                                 //
//                                                                            //
// with                                                                       //
//                  R0   R            pi*R                                    //
//   Sig(R) = exp(- -- - -- - eps cos(-----))                                 //
//                  R    Rd            Rd                                     //
//                                                                            //
// and                                                                        //
//                                                                            //
//   h(z)   = delta(z)                   for  d=0                             //
//            (1/2 d)  * exp(-|z/d|)     for  d>0                             //
//            (1/4|d|) * sech^2(|z/2d|)  for  d<0                             //
//                                                                            //
// The potential part returned by operator() amounts to                       //
//                                                                            //
//                                                                            //
//   Phi(r,z) = Sig(r) * H(z)                                                 //
//                                                                            //
// where r = sqrt(R^2+z^2) is the spherical polar radius and H(z) is the      //
// second integral of h(z).                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void DiskAnsatz::setup(const DiskPar& p)
{
  S0 = abs(p(0));
  Rd = abs(p(1));
  zd = abs(p(2));
  R0 = abs(p(3));
  eps= p(4);
  if(zd==0.)  thin       = 1;
  else    thin       = 0;
  if(R0==0.)  hollow     = 0;
  else    hollow     = 1;
  if(p(2)<0.) isothermal = 1;
  else    isothermal = 0;
  zdoRd = zd/Rd;
  R0oRd = R0/Rd;
  Rd2   = Rd*Rd;
  fac   = TPi*Units::Grav*S0;
}

double DiskAnsatz::SurfaceDensity(const double R) const
{
  if(hollow && R==0.) return 0.;
  register double y=R/Rd;
  if(eps==0.)         return S0*exp(-R0/R-y);
  return S0*exp(-R0/R-y+eps*cos(Pi*y));
}

inline double DiskAnsatz::mass_integrand(const double y) const
{
  if(y<=0. || y>=1.) return 0.;
  register double y1=1.-y, x=y/y1;    
  if(eps) return exp(-R0oRd/x-x+eps*cos(Pi*x))*x/(y1*y1);
  return exp(-R0oRd/x-x)*x/(y1*y1);
}

double DiskAnsatz::mass(const double R) const
{
  register double F=TPi*S0*Rd2;
  if(R<=0.) {         // give total mass
    if(eps)
      return F*qbulir(this,&DiskAnsatz::mass_integrand,0.,1.,1.e-6);
    if(hollow) return FPi*S0*R0*Rd*Kn(2,2.*sqrt(R0/Rd));
    return F;
  }                   // give mass inside R (but integrate z from -oo to oo)
  return F*qbulir(this,&DiskAnsatz::mass_integrand,0.,R/(Rd+R),1.e-6);
}

double DiskAnsatz::Density(const double R, const double z) const
{
  register double rhR;
  if(eps)          rhR= (hollow && R==0.)? 0. : exp(-R0/R-R/Rd+eps*cos(Pi*R/Rd));
  else if(hollow)  rhR= (R==0.)? 0. : exp(-R0/R-R/Rd);
  else             rhR= exp(-R/Rd);
  if(thin)                            // vertically thin disk: return Sigma
    return S0*rhR;
  else if(isothermal) {               // vertically isothermal disk
    register double x  =abs(z/zd),
      ex =exp(-x),
      ex1=1.+ex;
    return S0*rhR * ex/(ex1*ex1*zd);
  }                                   // vertically exponential disk
  register double x  =abs(z/zd),
    ex =exp(-x);
  return S0*rhR * 0.5*ex/zd;
}

double DiskAnsatz::Residual(const double r, const double st, const double ct)
  const
// gives aimed Laplace(Phi_multipole)
{
  if(ct==0. || S0==0.) return 0;
  register double R=r*st, z=r*ct, g,gp,gpp, F,f,fp,fpp;
  // deal with the vertical part
  if(thin) {
    g   = abs(z);
    gp  = sign(z);
    gpp = 0.;
  } else if(isothermal) {
    register double x,sh1;
    x   = abs(z/zd);
    gpp = exp(-x);
    sh1 = 1.+gpp;
    g   = 2*zd*(0.5*x+log(0.5*sh1));
    gp  = sign(z)*(1.-gpp)/sh1;
    gpp/= 0.5*sh1*sh1*zd;
  } else {
    register double x;
    x   = abs(z/zd);
    gpp = exp(-x);
    g   = zd*(gpp-1+x);
    gp  = sign(z)*(1.-gpp);
    gpp/= zd;
  }
  // deal with the radial part
  if(hollow && r==0.) F=f=fp=fpp=0.;
  else if(eps) {
    if(hollow) {
      register double rq=r*r,cr=cos(Pi*r/Rd),sr=sin(Pi*r/Rd);
      F   = (R==0)? 0. : exp(-R0/R-R/Rd+eps*cos(Pi*R/Rd));
      f   = exp(-R0/r-r/Rd+eps*cr);
      fp  = R0/rq-(1.+eps*sr)/Rd;
      fpp = (fp*fp-2.*R0/(rq*r)-eps*cr/Rd2)*f;
      fp *= f;
    } else {
      register double cr=cos(Pi*r/Rd),sr=sin(Pi*r/Rd);
      F   = exp(-R/Rd+eps*cos(Pi*R/Rd));
      f   = exp(-r/Rd+eps*cr);
      fp  = -(1.+eps*sr)/Rd;
      fpp = (fp*fp-eps*cr/Rd2)*f;
      fp *= f;
    }
  } else {
    if(hollow) {
      register double rq=r*r;
      F   = (R==0)? 0. : exp(-R0/R-R/Rd);
      f   = exp(-R0/r-r/Rd);
      fp  = R0/rq-1./Rd;
      fpp = (fp*fp-2.*R0/(rq*r))*f;
      fp *= f;
    } else {
      F   = exp(-R/Rd);
      f   = exp(-r/Rd);
      fp  =-f/Rd;
      fpp = f/Rd2;
    }
  }

  return fac * ((F-f)*gpp - 2*fp*(g+z*gp)/r - fpp*g);
}

double DiskAnsatz::operator() (const double R, const double z, const double r,
                               double* dP) const
{
  if(S0==0.) {
    if(dP) dP[0]=dP[1]=0.;
    return 0.;
  }
  register double g,f;
  if(dP) {
    register double gp,fp;
    if(thin) {
      g  = abs(z);
      gp = sign(z);
    } else if(isothermal) {
      register double x,ex,sh1;
      x  = abs(z/zd);
      ex = exp(-x);
      sh1= 1.+ex;
      g  = 2*zd*(0.5*x+log(0.5*sh1));
      gp = sign(z)*(1.-ex)/sh1;
    } else {
      register double x,ex;
      x  = abs(z/zd);
      ex = exp(-x);
      g  = zd*(ex-1+x);
      gp = sign(z)*(1.-ex);
    }

    if(hollow && r==0.) f=fp=0.;
    else if(eps) {
      if(hollow) {
	register double rq=r*r,cr=cos(Pi*r/Rd),sr=sin(Pi*r/Rd);
	f   = exp(-R0/r-r/Rd+eps*cr);
	fp  = (R0/rq-(1.+eps*sr)/Rd)*f;
      } else {
	register double cr=cos(Pi*r/Rd),sr=sin(Pi*r/Rd);
	f   = exp(-r/Rd+eps*cr);
	fp  = -(1.+eps*sr)*f/Rd;
      }
    } else {
      if(hollow) {
	register double rq=r*r;
	f   = exp(-R0/r-r/Rd);
	fp  = (R0/rq-1./Rd)*f;
      } else {
	f  = exp(-r/Rd);
	fp =-f/Rd;
      }
    }
    dP[0] = fac *  R/r * fp * g;
    dP[1] = fac * (z/r * fp * g + f * gp);
  } else {
    if(thin) {
      g  =abs(z);
    } else if(isothermal) {
      register double x,ex,sh1;
      x  = abs(z/zd);
      ex = exp(-x);
      sh1= 1.+ex;
      g  = 2*zd*(0.5*x+log(0.5*sh1));
    } else {
      register double x,ex;
      x  = abs(z/zd);
      ex = exp(-x);
      g  = zd*(ex-1+x);
    }

    if(hollow && r==0.) f=0.;
    else if(eps) {
      if(hollow) f = exp(-R0/r-r/Rd+eps*cos(Pi*r/Rd));
      else       f = exp(-r/Rd+eps*cos(Pi*r/Rd));
    } else {
      if(hollow) f = exp(-R0/r-r/Rd);
      else       f = exp(-r/Rd);
    }
  }

  return fac*f*g;
}

double DiskAnsatz::Laplace(const double R, const double z) const
{
  if(S0==0.) return 0;
  register double r=hypot(R,z), g,gp,gpp, F,f,fp,fpp;
  // deal with the vertical part
  if(thin) {
    g  =abs(z);
    gp =sign(z);
    gpp=0.;
  } else if(isothermal) {
    register double x,sh1;
    x   = abs(z/zd);
    gpp = exp(-x);
    sh1 = 1.+gpp;
    g   = 2*zd*(0.5*x+log(0.5*sh1));
    gp  = sign(z)*(1.-gpp)/sh1;
    gpp/= 0.5*sh1*sh1*zd;
  } else {
    register double x;
    x   = abs(z/zd);
    gpp = exp(-x);
    g   = zd*(gpp-1+x);
    gp  = sign(z)*(1.-gpp);
    gpp/= zd;
  }
  // deal with the radial part
  if(hollow && r==0.) F=f=fp=fpp=0.;
  else if(eps) {
    if(hollow) {
      register double rq=r*r,cr=cos(Pi*r/Rd),sr=sin(Pi*r/Rd);
      F   = (R==0)? 0. : exp(-R0/R-R/Rd+eps*cos(Pi*R/Rd));
      f   = exp(-R0/r-r/Rd+eps*cr);
      fp  = R0/rq-(1.+eps*sr)/Rd;
      fpp = (fp*fp-2.*R0/(rq*r)-eps*cr/Rd2)*f;
      fp *= f;
    } else {
      register double cr=cos(Pi*r/Rd),sr=sin(Pi*r/Rd);
      F   = exp(-R/Rd+eps*cos(Pi*R/Rd));
      f   = exp(-r/Rd+eps*cr);
      fp  = -(1.+eps*sr)/Rd;
      fpp = (fp*fp-eps*cr/Rd2)*f;
      fp *= f;
    }
  } else {
    if(hollow) {
      register double rq=r*r;
      F   = (R==0)? 0. : exp(-R0/R-R/Rd);
      f   = exp(-R0/r-r/Rd);
      fp  = R0/rq-1./Rd;
      fpp = (fp*fp-2.*R0/(rq*r))*f;
      fp *= f;
    } else {
      F   = exp(-R/Rd);
      f   = exp(-r/Rd);
      fp  =-f/Rd;
      fpp = f/Rd2;
    }
  }
  return fac * (f*gpp + 2*fp*(g+z*gp)/r + fpp*g);
}

Frequencies DiskAnsatz::kapnuom(const double R) const
// returns dPhi/dR, d^2Phi/dR^2, d^2Phi/dz^2 at z=0
{
  if(S0==0.) return Frequencies(0.);
  register double er, gpp;
  if(hollow) er = (R==0)? 0. : exp(-R0/R-R/Rd+eps*cos(Pi*R/Rd));
  else   er = exp(-R/Rd+eps*cos(Pi*R/Rd));
  if(thin) {                          // vertically thin disk
    std::cerr<<"Warning: KapNuOm(Phi) involves delta(z) at z=0\n";
    gpp=0.;
  } else if(isothermal)               // vertically isothermal disk
    gpp=0.5/zd;
  else                                // vertically exponential disk
    gpp=1./zd;
  Frequencies om;
  om[0] = 0.;
  om[1] = fac*er*gpp;
  om[2] = 0.;
  return om;
}

////////////////////////////////////////////////////////////////////////////////
// class Disks
////////////////////////////////////////////////////////////////////////////////

void Disks::reset(const int N, const DiskPar* p)
{
  delete[] D;
  nd = N;
  register int i;
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  for(i=0; i<nd; i++) (D+i)->setup(p[i]);
}

Disks::Disks(std::istream& from)
{
  if(!from) {
    std::cerr<<" Trying to construct Disks from a closed std::istream\n";
    std::exit(1);
  }
  DiskPar P;
  from >> nd;
  SwallowRestofLine(from);
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  register DiskAnsatz *p=D;
  for(; p<Dup; p++) {
    from >> P;
    SwallowRestofLine(from);
    p->setup(P);
  }
  //  nemo_dprintf(4,"Disks: read %d parameters\n",nd);
}

Disks::Disks(const Disks& DS) : nd(DS.nd)
{
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  register DiskAnsatz *p=D;
  for(register int i=0; p<Dup; p++,i++)
    p->setup(DS.Parameter(i));
}

Disks::Disks(const int N, const DiskPar* p) : nd(N)
{
  register int i;
  D   = new DiskAnsatz[nd];
  Dup = D+nd;
  for(i=0; i<nd; i++) (D+i)->setup(p[i]);
}

////////////////////////////////////////////////////////////////////////////////
// class SpheroidDensity
////////////////////////////////////////////////////////////////////////////////

void SpheroidDensity::setup(const SphrPar& d)
{
  rh0 = d(0);
  q   = d(1);
  gam = d(2);
  bet = d(3);
  r0  = d(4);
  rcut= d(5);

  if(rcut<=0.) rci = 0.;
  else     rci = 1./rcut;
  beg = bet-gam;
  qi  = 1./q;
  r0i = 1./r0;
}

double SpheroidDensity::Density(const double R, const double z) const
{
  register double m = hypot(R,z*qi), m0=m*r0i, rho=rh0;
  if(gam==0.5)   rho /= sqrt(m0);
  else if(gam==1.)    rho /= m0;
  else if(gam==2.)    rho /= m0*m0;
  else if(gam!=0.)    rho /= pow(m0,gam);
  m0 += 1;
  if(beg==1.)    rho /= m0;
  else if(beg==2.)    rho /= m0*m0;
  else if(beg==3.)    rho /= m0*m0*m0;
  else                rho /= pow(m0,beg);
  if(rci)             rho *= exp(-square(m*rci));
  return rho;
}

double SpheroidDensity::mass_integrand(const double y) const
{
  if(rci) {
    register double y1=1.-y, m=r0*y/y1;
    return pow(y,2.-gam) * pow(y1,bet-4.) * exp(-square(m*rci));
  }
  return pow(y,2-gam) * pow(1.-y,bet-4.);
}

double SpheroidDensity::mass(const double m) const
{
  return FPi*q*rh0*cube(r0)*
    qbulir(this,&SpheroidDensity::mass_integrand,0.,m/(m+r0),1.e-6);
}

////////////////////////////////////////////////////////////////////////////////
// class Spheroids
////////////////////////////////////////////////////////////////////////////////

void Spheroids::reset(const int N, const SphrPar* p)
{
  delete[] S;
  ns = N;
  register int i;
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  for(i=0; i<ns; i++) (S+i)->setup(p[i]);
}

Spheroids::Spheroids(std::istream& from)
{
  if(!from) {
    std::cerr<<" Trying to construct Spheroids from a closed istream\n";
    std::exit(1);
  }
  SphrPar P;
  from >> ns;
  SwallowRestofLine(from);
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  SpheroidDensity *p=S;
  for(; p<Sup; p++) {
    from >> P;
    SwallowRestofLine(from);
    p->setup(P);
  }
  //nemo_dprintf(4,"Spheroids: read %d parameters\n",ns);
}

Spheroids::Spheroids(const Spheroids& SP) : ns(SP.ns)
{
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  SpheroidDensity *p=S;
  for(register int i=0; p<Sup; p++,i++) 
    p->setup(SP.Parameter(i));
}

Spheroids::Spheroids(const int N, const SphrPar* p) : ns(N)
{
  register int i;
  S   = new SpheroidDensity[ns];
  Sup = S+ns;
  for(i=0; i<ns; i++) (S+i)->setup(p[i]);
}

double Spheroids::beta() const
{
  register double b=1.e3;
  for(register SpheroidDensity *p=S; p<Sup; p++)
    b = min(b, p->outer_power());
  return (b==1.e3)? -1 : b;
}

double Spheroids::gamma() const
{
  register double g=0.;
  for(register SpheroidDensity *p=S; p<Sup; p++)
    g = max(g, p->inner_power());
  return g;
}

////////////////////////////////////////////////////////////////////////////////
// class Multipole
////////////////////////////////////////////////////////////////////////////////
const int N =LMAX/2+1;   // number of multipoles
const int N2=3*N/2;	   // number of grid point for cos[theta] in [0,1]
const int N4=5*N/2; 	   // number of points used to integrate over cos[theta]

typedef Vector<double,N> DBN;

void Multipole::AllocArrays()
{
  if(LR) {
    lLc = new double[K[0]];
    d2R = new double[K[0]];
    d2L = new double[K[0]];
  }
  logr = new double[K[0]];
  X[0] = logr;
  X[1] = new double[K[1]];
  Alloc2D(Y[0],K); Alloc2D(Y[1],K); Alloc2D(Y[2],K); // See FreeMemory.h
  Alloc2D(Z[0],K); Alloc2D(Z[1],K); Alloc2D(Z[2],K); Alloc2D(Z[3],K);
}

Multipole::Multipole(const int Kk,
                     const double ri,
                     const double ra,
                     const double g,
                     const double b,
                     PotResidual const *PR,
                     const int    lr)
{
  //nemo_dprintf(4,"Multipole::Multipole() ... \n");
  LR   = lr;
  K[0] = Kk;
  K[1] = N2;
  AllocArrays();
  setup(ri,ra,g,b,PR);
  //nemo_dprintf(4," done Multipole::Multipole()\n");
}

void Multipole::reset(const double ri,
                      const double ra,
                      const double g,
                      const double b,
                      PotResidual const *PR,
                      const int lr)
{
  if(lr && !LR) {    // so it doesn't crash
    lLc = new double[K[0]];
    d2R = new double[K[0]];
    d2L = new double[K[0]];
  }
  if(!lr && LR) {    // to avoid memory leaks
    delete[] lLc;
    delete[] d2R;
    delete[] d2L;
  }
  LR = lr;
  setup(ri,ra,g,b,PR);
}

void Multipole::setup(const double ri, const double ra,
                      const double g, const double b,
                      PotResidual const *PR)
{
  Rmin = ri;
  Rmax = ra;
  gamma= g; 
  beta = b;
  lRmin= log(Rmin);
  lRmax= log(Rmax); 
  g2   = 2.-gamma;

  const    DBN    Zero=DBN(0.);
  const    double half=0.5, three=3., sixth=1./6.,
    dlr =(lRmax-lRmin)/double(K[0]-1);
  register int    i,l,k,ll,lli1;
  register double dx,dx2,xl_ll,xh_ll,risq,ril2,dP;
  DBN    A[4],P2l,dP2l;

  register DBN    EX;
  //
  // 0  check for inconsistencies in input
  //
  //nemo_dprintf(5,"Multipole::setup(): 0\n");
  if(beta>0. && beta<3.) {
    std::cerr<<" Warning: beta= "<<beta
	     <<" unsuitable for Multipole expansion;"
	     <<" we'll take beta=3.2\n";
    beta=3.2;
  }
  //
  // 1  compute expansion of the density
  //
  //nemo_dprintf(5,"Multipole::setup(): 1\n");
  double
    *ct   = new double[N4],
    *st   = new double[N4],
    *wi   = new double[N4],
    *r    = new double[K[0]];
  DBN
    *W    = new DBN   [N4],
    *rhol = new DBN   [K[0]],
    *rhl2 = new DBN   [K[0]];
  //
  // 1.1 set points and weights for integration over cos(theta)
  //
  //nemo_dprintf(6,"Multipole::setup(): 1.1\n");
  GaussLegendre(ct,wi,N4);
  //nemo_dprintf(8,"Multipole::setup(): N4=%d\n",N4);
  for(i=0; i<N4; i++) {
    ct[i] = 0.5 * (ct[i]+1.);
    st[i] = sqrt(1.-ct[i]*ct[i]);
    wi[i] = 0.5 * wi[i];
    LegendrePeven(W[i],ct[i]);
    W[i] *= wi[i];
  }
  //
  // 1.2 integrate over cos(theta)
  //
  //nemo_dprintf(6,"Multipole::setup(): 1.2\n");
  for(k=0; k<K[0]; k++) {
    logr[k] = k<K[0]-1? lRmin+dlr*k : lRmax;                 // v0.7
    r[k]    = exp(logr[k]);
    rhol[k] = 0.;
    for(i=0; i<N4; i++)
      rhol[k] += W[i] * PR->Residual(r[k],st[i],ct[i]);
  }
  delete[] ct;
  delete[] st;
  delete[] wi;
  delete[] W;
  //
  // 1.3 establish spline in r needed for integration
  //
  //nemo_dprintf(6,"Multipole::setup(): 1.3\n");
  spline(r,rhol,K[0],(-gamma/r[0])*rhol[0],Zero,rhl2,0,1);
  //
  // 2. compute potential's expansion
  //
  //nemo_dprintf(5,"Multipole::setup(): 2\n");
  DBN
    *P1   = new DBN[K[0]],
    *P2   = new DBN[K[0]],
    *Phil = new DBN[K[0]],
    *dPhl = new DBN[K[0]];
  //
  // 2.1 set P1[k][l] r[k]^(-1-2l) = Int[rho_2l(x,l) x^(2l+2), {x,0,r[k]}]
  //
  //     for r < Rmin we take  rho_2l proportional r^-gamma
  //
  //nemo_dprintf(6,"Multipole::setup(): 2.1\n");
  risq  = Rmin*Rmin;
  for(l=0; l<N; l++) {
    P1[0][l] = rhol[0][l] * risq / double(2*l+3-gamma);
    EX[l]    = exp(-(1+2*l)*dlr);
  }
  for(k=0; k<K[0]-1; k++) {
    dx   = r[k+1]-r[k];
    dx2  = dx*dx;
    A[0] = r[k+1]*rhol[k] - r[k]*rhol[k+1] + sixth*r[k]*r[k+1] *
      ( (r[k+1]+dx)*rhl2[k] - (r[k]-dx)*rhl2[k+1] );
    A[1] = rhol[k+1]-rhol[k]
      + sixth * ( (dx2-three*r[k+1]*r[k+1]) * rhl2[k]
		  -(dx2-three*r[k]*r[k])     * rhl2[k+1] );
    A[2] = half  * (r[k+1]*rhl2[k] - r[k]*rhl2[k+1]);
    A[3] = sixth * (rhl2[k+1]-rhl2[k]);
    for(l=0,ll=2; l<N; l++,ll+=2) {
      xl_ll = r[k]*EX(l);
      xh_ll = r[k+1];
      for(i=0,lli1=ll+1,dP=0.; i<4; i++,lli1++) {
	xl_ll*= r[k];
	xh_ll*= r[k+1];
	dP   += A[i](l) * (xh_ll - xl_ll) / lli1;
      }
      P1[k+1][l] = EX(l) * P1[k](l) + dP / dx;
    }
  }
  //
  // 2.2 set P2[k][l] = r[k]^(2l) Int[rho_2l(x,l) x^(1-2l), {x,r[k],Infinity}]
  //
  //     for r > Rmax we take  rho_2l proportional r^-beta if beta>0
  //                                  = 0                  if beta<=0
  //
  //nemo_dprintf(6,"Multipole::setup(): 2.2\n");
  if(beta>0.) {
    risq  = Rmax*Rmax;
    for(l=0; l<N; l++) {
      P2[K[0]-1][l] = rhol[K[0]-1][l] * risq / double(beta+2*l-2);
      EX[l] = exp(-2*l*dlr);
    }
  } else {
    P2[K[0]-1] = 0.;
    for(l=0; l<N; l++)
      EX[l] = exp(-2*l*dlr);
  }
  for(k=K[0]-2; k>=0; k--) {
    risq = r[k]*r[k];
    dx   = r[k+1]-r[k];
    dx2  = dx*dx;
    A[0] = r[k+1]*rhol[k] - r[k]*rhol[k+1] + sixth*r[k]*r[k+1] *
      ( (r[k+1]+dx)*rhl2[k] - (r[k]-dx)*rhl2[k+1] );
    A[1] = rhol[k+1]-rhol[k]
      + sixth * ( (dx2-three*r[k+1]*r[k+1]) * rhl2[k]
		  -(dx2-three*r[k]*r[k])     * rhl2[k+1] );
    A[2] = half  * (r[k+1]*rhl2[k] - r[k]*rhl2[k+1]);
    A[3] = sixth * (rhl2[k+1]-rhl2[k]);
    for(l=0,ll=1,ril2=1.; l<N; l++,ll-=2,ril2*=risq) {
      xl_ll = r[k];
      xh_ll = r[k+1]*EX(l);
      for(i=0,lli1=ll+1,dP=0.; i<4; i++,lli1++) {
	xl_ll *= r[k];
	xh_ll *= r[k+1];
	if(lli1) dP += A[i](l) * (xh_ll - xl_ll) / lli1;
	else     dP += A[i](l) * ril2 * dlr;
      }
      P2[k][l] = EX(l) * P2[k+1](l) + dP / dx;
    }
  }
  //
  // 2.3 put together the Phi_2l(r) and dPhi_2l(r)/dlog[r]
  //
  //nemo_dprintf(6,"Multipole::setup(): 2.3\n");
  for(k=0; k<K[0]; k++)
    for(l=ll=0; l<N; l++,ll+=2) {
      Phil[k][l] =-P1[k](l) - P2[k](l);                   // Phi_2l
      dPhl[k][l] = (ll+1)*P1[k](l) - ll*P2[k](l);         // dPhi_2l/dlogr
    }
  if(gamma<2)
    Phi0 = Phil[0](0) - dPhl[0](0) / g2;
  delete[] r;
  delete[] rhol;
  delete[] rhl2;
  delete[] P1;
  delete[] P2;
  //
  // 3. establish L_circ(R) on the logarithmic grid
  //
  //nemo_dprintf(5,"Multipole::setup(): 3\n");
  if(LR) {
    tg3 = 2./(3.-gamma);
    g3h = 0.5*(3.-gamma);
    LegendrePeven(P2l,0.);
    for(k=0; k<K[0]; k++)
      lLc[k] = 0.5 * ( 2*logr[k] + log(dPhl[k]*P2l) );
    spline(lLc,logr,K[0],tg3,2.,d2R,0,0);
    spline(logr,lLc,K[0],g3h,.5,d2L,0,0);
    lzmin = lLc[0];
    lzmax = lLc[K[0]-1];
  }
  //
  // 4.  Put potential and its derivatives on a 2D grid in log[r] & cos[theta]
  //
  //nemo_dprintf(5,"Multipole::setup(): 4\n");
  //
  // 4.1 set linear grid in theta
  //
  //nemo_dprintf(6,"Multipole::setup(): 4.1\n");
  for(i=0; i<N2; i++) 
    X[1][i] = double(i) / double(N2-1);
  //
  // 4.2 set dPhi/dlogr & dPhi/dcos[theta] 
  //
  //nemo_dprintf(6,"Multipole::setup(): 4.2\n");
  for(i=0; i<N2; i++) {
    dLegendrePeven(P2l,dP2l,X[1][i]);
    for(k=0; k<K[0]; k++) {
      Y[0][k][i] = Phil[k] * P2l;			// Phi
      Y[1][k][i] = dPhl[k] * P2l;			// d Phi / d logR
      Y[2][k][i] = Phil[k] * dP2l;		// d Phi / d cos(theta)
    }
  }
  delete[] Phil;
  delete[] dPhl;
  //
  // 4.3 establish 2D Pspline of Phi in log[r] & cos[theta]
  //
  //nemo_dprintf(6,"Multipole::setup(): 4.3\n");
  Pspline2D(X,Y,K,Z);
  //nemo_dprintf(5,"Multipole::setup(): done\n");
}

Multipole::~Multipole()
{
  if(LR) {
    delete[] lLc;
    delete[] d2R;
    delete[] d2L;
  }
  delete[] X[0]; delete[] X[1];
  Free2D(Y[0]); Free2D(Y[1]); Free2D(Y[2]);
  Free2D(Z[0]); Free2D(Z[1]); Free2D(Z[2]); Free2D(Z[3]);
}

double Multipole::operator() (const double r, const double ct, const double st,
                              double* dP) const
{
  double Xi[2];
  register double lr=log(r), Phi;
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = abs(ct);
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP);
  if(dP) dP[1]*= sign(ct);
  if(lr < lRmin) {
    if(g2>0.) {
      Phi = (Phi-Phi0)*exp(g2*(lr-Xi[0]));
      if(dP) dP[0] = g2*Phi;
      Phi+= Phi0;
    } else if(g2==0.) {
      if(dP) dP[0] = Phi/lRmin;
      Phi*= lr/lRmin;
    } else {
      Phi*= exp(g2*(lr-Xi[0]));
      if(dP) dP[0] = g2*Phi;
    }
  } else if(lr > lRmax) {
    Phi *= Rmax/r;
    if(dP) dP[0] =-Phi;
  }
  if(dP) {
    register double temp;
    dP[0]/= r;
    dP[1]*=-st/r;
    temp  = ct*dP[0] - st*dP[1];
    dP[0] = st*dP[0] + ct*dP[1];
    dP[1] = temp;
  }
  return Phi;
}

double Multipole::vcsquare(const double R, double &dvcqdR) const
{
  const int n2[2]={2,2};
  double Xi[2], dP[2], **d2P;
  register double lr=log(R), Phi;
  Alloc2D(d2P,n2);
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = 0.;
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP,d2P);
  if(lr < lRmin) {
    if(g2>0.) {
      dP[0]     = g2*(Phi-Phi0)*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    } else if(g2==0.) {
      dP[0]     = Phi/lRmin;
      d2P[0][0] = 0.;
    } else {
      dP[0]     = g2*Phi*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    }
  } else if(lr > lRmax) {
    dP[0]     =-Phi*Rmax/R;
    d2P[0][0] =-dP[0];
  }
  dvcqdR = d2P[0][0] / R;     // dvc^2 / dR = (1/R) d^2 Phi/ d(lnR)^2
  Free2D(d2P);
  return dP[0];               // vc^2       = dPhi/dlnR
}

double Multipole::RfromLc (const double L, double* dR) const
{
    if(!LR || L<=0.) {
        if(!LR) cerr<<" Multipole.RfromLc() was not initialized\n";
        if(dR) *dR=0.;
	return 0.;
    }
    register double lL=log(L);
    if(dR) {
        register double R;
        if(lL<lzmin) {
	    *dR = tg3;
	    R   = Rmin * exp(tg3*(lL-lzmin));
        } else if(lL>lzmax) {
	    *dR = 2.;
	    R   = Rmax * exp(2.*(lL-lzmax));
        } else
            R = exp(splev(lLc,logr,d2R,K[0],lL,dR));
        *dR *= R/L;
	return R;
    }
    if(lL<lzmin)      return Rmin * exp(tg3*(lL-lzmin));
    else if(lL>lzmax) return Rmax * exp(2.*(lL-lzmax));
                      return exp(splev(lLc,logr,d2R,K[0],lL));
}

double Multipole::LfromRc (const double R, double* dL) const
{
    if(!LR || R<=0.) {
        if(!LR) cerr<<" Multipole.LfromRc() was not initialized\n";
        if(dL) *dL=0.;
	return 0.;
    }
    register double lR=log(R);
    if(dL) {
        register double L;
        if(lR<lRmin) {
	    *dL = g3h;
	    L   = exp(lzmin + g3h*(lR-lRmin));
        } else if(lR>lRmax) {
	    *dL = 0.5;
	    L   = exp(lzmax + 0.5*(lR-lRmax));
        } else
            L = exp(splev(logr,lLc,d2L,K[0],lR,dL));
        *dL *= L/R;
	return L;
    }
    if(lR<lRmin)      return exp(lzmin + g3h*(lR-lRmin));
    else if(lR>lRmax) return exp(lzmax + 0.5*(lR-lRmax));
            	      return exp(splev(logr,lLc,d2L,K[0],lR));
}

double Multipole::Laplace(const double r, const double ct) const
{
  const    int    m[2]={2,2};
  register double lr=log(r), Phi, Lap;
  double Xi[2], *dP, **d2P;
  Alloc1D(dP,2);
  Alloc2D(d2P,m);
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = abs(ct);
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP,d2P);
  dP[1]*=sign(ct);
  if(lr < lRmin) {
    if(g2>0.) {
      dP[0]     = g2*(Phi-Phi0)*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    } else if(g2==0.) {
      dP[0]     = Phi/lRmin;
      d2P[0][0] = 0.;
    } else {
      dP[0]     = g2*Phi*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    }
  } else if(lr > lRmax) {
    dP[0]     =-Phi*Rmax/r;
    d2P[0][0] =-dP[0];
  }
  Lap = ( dP[0]+d2P[0][0] + (1.-ct*ct)*d2P[1][1]-2.*ct*dP[1] ) / (r*r);
  Free1D(dP);
  Free2D(d2P);
  return Lap;
}

Frequencies Multipole::kapnuom(const double R) const
// returns dPhi/dR, d^2Phi/dR^2/ d^2Phi/dz^2
{
  const    int    m[2]={2,2};
  register double lr=log(R), Rq=R*R, Phi;
  double Xi[2], *dP, **d2P;
  Alloc1D(dP,2);
  Alloc2D(d2P,m);
  Xi[0] = min(lRmax,max(lRmin,lr));
  Xi[1] = 0.;
  Phi   = Psplev2D(X,Y,Z,K,Xi,dP,d2P);
  if(lr < lRmin) {
    if(g2>0.) {
      dP[0]     = g2*(Phi-Phi0)*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    } else if(g2==0.) {
      dP[0]     = Phi/lRmin;
      d2P[0][0] = 0.;
    } else {
      dP[0]     = g2*Phi*exp(g2*(lr-Xi[0]));
      d2P[0][0] = g2*dP[0];
    }
  } else if(lr > lRmax) {
    dP[0]     =-Phi*Rmax/R;
    d2P[0][0] =-dP[0];
  }
  Frequencies om;
  om[0] = (d2P[0][0]-dP[0]) / Rq;
  om[1] = (dP[0] + d2P[1][1]) / Rq;
  om[2] = dP[0] / R;
  Free1D(dP);
  Free2D(d2P);
  return om; 
}

////////////////////////////////////////////////////////////////////////////////
// class GalaxyPotential
////////////////////////////////////////////////////////////////////////////////

double GalaxyPotential::operator() (const double R, const double z) const
{
  register double r  =hypot(R,z), pot=(r)? M(r,z/r,R/r) : M(0,0,0);
  for(register DiskAnsatz *p=D; p<Dup; p++) pot+= (*p)(R,z,r);
  return pot;
}

double GalaxyPotential::operator() (const double R, const double z,
                                    double& dR, double &dz) const
{
  double d[2];
  register double r  =hypot(R,z), pot=(r)? M(r,z/r,R/r,d) : M(0,0,0); 
  dR = (r)? d[0] : 0.; dz = (r)? d[1] : 0.;
  for(register DiskAnsatz *p=D; p<Dup; p++) {
    pot += (*p)(R,z,r,d);
    dR  += (r)? d[0] : 0.;
    dz  += (r)? d[1] : 0.;
  }
  return pot;
}

void GalaxyPotential::OortConstants(const double R, double &A, double &B) const
{
  double vc, dvc;
  vc  = sqrt(M.vcsquare(R,dvc));
  dvc/= 2.*vc;
  A   = 0.5 * (vc/R - dvc);
  B   =-0.5 * (vc/R + dvc);
}

double GalaxyPotential::Laplace(const double R, const double z) const
{
  register double r=hypot(R,z), L=M.Laplace(r,z/r);
  register DiskAnsatz *p=D;
  for(; p<Dup; p++) L += p->Laplace(R,z);
  return L;
}

Frequencies GalaxyPotential::KapNuOm  (const double R) const
{
  Frequencies Om = M.kapnuom(R);
  for(register DiskAnsatz *p=D; p<Dup; p++) Om += p->kapnuom(R);
  Om[2]/= R;
  Om[0]+= 3*Om(2);
  Om.apply(&sqrt);
  return Om;                  // omega, kappa, nu
}

#endif  // #ifndef GalPot_cc
