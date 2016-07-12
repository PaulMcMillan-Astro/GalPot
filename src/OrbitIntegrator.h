/***************************************************************************//**
\file  OrbitIntegrator.h
\brief Classes and functions used for integrating orbits in gravitational potentials

*                                                                              *
* C++ code written by Walter Dehnen, 1994-96,                                  *
*                     Paul McMillan, 2007-                                     *
* e-mail: paul@astro.lu.se                                                     *
* github: https://github.com/PaulMcMillan-Astro/Torus                          *
*                                                                              *
*******************************************************************************/

#ifndef _TorusOrb_
#define _TorusOrb_ 1

#include "Vector.h"
#include "Potential.h"
#include "Numerics.h"


/** \brief Integrates orbits in gravitational potentials. */

class OrbitIntegratorStep {
private:
    Potential*  Pot;                       // pointer to the potential
    Vector <double,6> QP;  // R,z, phi, vR,vz, vphi

    double      Jphi, E, tol, dtm,dtmin;   // energy, tolerance, max time step, min time step
    void        RungeKutta(const double dt);  // 4th order Runge-Kutta integrator
    Vector<double,6> RK3Dint(const Vector<double,6> Y);
public:
    OrbitIntegratorStep           () {;}
    OrbitIntegratorStep           (const Vector <double,6>, Potential*, const double=0.);
    void setup                    (const Vector <double,6>, Potential*, const double=0.);
    void setup                    (const Vector <double,6>);
   ~OrbitIntegratorStep           () {;}
    void       set_tolerance(const double dE=0.) {tol= (dE)? fabs(dE) : 1.e-9;}
    void       set_maxstep  (const double dt=0.) {dtm= (dt)? fabs(dt) : 1.;}
    void       stepRK_by    (double&, const double=1.4);
    
    Vector<double,6> XV        () const {
    Vector<double,6> xv=QP; xv[5] = xv[0]*xv[5]; return xv; }
    
    double     Energy()     {return E;}
    double     AngularMomentum()     {return Jphi;}
    double     maxstep()    {return dtm;}
};

// Class that integrates orbit for a time T, determining min, max and mean R, max z, (SOS in the future)?

class OrbitIntegratorWithStats {
  Potential *Pot;
  OrbitIntegratorStep Stepper;
  bool setupDone;
  double Tmax;
  int runGeneric(const string type,Vector <double,6> *output,double *tout, int N);
public:
  OrbitIntegratorWithStats() {;}
  OrbitIntegratorWithStats(Vector<double,6> StartPoint, Potential *PotIn,
			   const double TmaxIn=13800.);
  ~OrbitIntegratorWithStats() {;}
  
  void setup(Vector<double,6> StartPoint, Potential *PotIn, const double TmaxIn);
  void setup(Vector<double,6> StartPoint);

  Vector<double,6> XV_ini;
  double Energy;
  double Lz;
  double MinR;
  double MaxR;
  double Minr;
  double Maxr;
  double Maxz;
  double GuidingRadius;
  double MeanR;
  double PseudoEccentricity;
  
  int run() { return runGeneric("NoOutput",NULL,NULL,1); }
  
  int runWithOutput(Vector<double,6>*out, int N)
    { return runGeneric("OutputNoTimes",out,NULL,N); }
  // Output N points between 0 and Ttot (roughly evenly spaced)
  
  int runWithOutputIncludingTime(Vector<double,6>*out, double *t, int N)
    { return runGeneric("OutputWithTimes",out,t,N); }
   // Output N points between 0 and Ttot, giving T in each case.
  
  int runWithOutputAtGivenTimes(Vector<double,6>*out, double*T,int N)
    { return runGeneric("OutputSetTimes",out,T,N); }
  // Output N points at times T
};


#endif
