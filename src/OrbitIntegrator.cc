/*******************************************************************************
*                                                                              *
*  OrbitIntegrator.cc                                                          *
*                                                                              *
* C++ code written by Paul McMillan, 2016                                      *
* e-mail: paul@astro.lu.se                                                     *
* github: https://github.com/PaulMcMillan-Astro                                *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <fstream>
#include "OrbitIntegrator.h"
using std::ofstream;

Vector<double,6> GiveVec6(const double &a0,const double &a1,const double &a2,
			  const double &a3,const double &a4,const double &a5) {
  Vector<double,6> out;
  out[0] = a0;
  out[1] = a1;
  out[2] = a2;
  out[3] = a3;
  out[4] = a4;
  out[5] = a5;
  return out;
}


typedef Vector<double,4>	DB4;
typedef Vector<double,6>	Vec6;
typedef Matrix<double,4,4>	DB44;

////////////////////////////////////////////////////////////////////////////////
OrbitIntegratorStep::OrbitIntegratorStep(const Vector<double,6> xv, Potential* Phi, const double dE)
{
  setup(xv,Phi,dE);
}

void OrbitIntegratorStep::setup(const Vector<double,6> xv, Potential* Phi, const double dE) 
 {
  Pot = Phi;
  set_tolerance(dE);
  setup(xv);
}

void OrbitIntegratorStep::setup(const Vector<double,6> xv) 
 {
  QP = xv;
  QP[5] = xv[5]/xv[0];
  Jphi = xv[5]*xv[0];
  Pot->set_Lz(Jphi);
  E   = Pot->eff(QP[0],QP[1]) + 0.5 * (QP[3]*QP[3] + QP[4]*QP[4]);
  // Could arrange something more automatic here
  set_maxstep();
  dtmin = 1e-10;
}



Vec6 OrbitIntegratorStep::RK3Dint(const Vector<double,6> Y)
{
  double aR,az,phidd;
  Pot->eff(Y(0),Y(1),aR,az);
  phidd = -2.*Y(3)*Y(5)/Y(0);
  return GiveVec6(Y(3),Y(4),Y(5),-aR,-az,phidd); // N.B. 3rd/6th component not 
                                             // velocity/acceleration (angular) 
}
 
void OrbitIntegratorStep::RungeKutta(const double dt)
{
  QP[5] = Jphi/(QP[0]*QP[0]); 
  Vec6 dY,Y0=QP,Y1=Y0;
  
  // Failure to reuse
  Y1+=(dY=dt*RK3Dint(Y0))/6.;
  Y1+=(dY=dt*RK3Dint(Y0+.5*dY))/3.;
  Y1+=(dY=dt*RK3Dint(Y0+.5*dY))/3.;
  Y1+=(dY=dt*RK3Dint(Y0+dY))/6.;
  //R=Y1(0); z=Y1(1); phi=Y1(2); pR=Y1(3); pz=Y1(4); phid=Y1(5);
  QP = Y1;
  E=Pot->eff(QP[0],QP[1])+0.5*(QP[3]*QP[3]+QP[4]*QP[4]);
  //phidd = -2.*Y1(3)*Y1(5)/Y1(0);
}

void OrbitIntegratorStep::stepRK_by(double& dt, const double f)
{
  //    register OrbitIntegratorStep next = *this;
  Vec6 oldQP = QP;
  double oldE = E;
  double fac=(f<=1.)? 1.4:f, FAC=pow(fac,5), dE;
  if(dt>dtm) dt = dtm;
  RungeKutta(dt);
  dE = fabs(oldE - E);
  while(dt*fac<dtm && FAC*dE < tol) {
    dt  *= fac;
    QP=oldQP; E = oldE;
    RungeKutta(dt);
    dE = fabs(E - oldE);
  }
  while(dE>tol && dt>dtmin*fac)  {
      dt  /= fac;
      QP=oldQP; E = oldE;
      RungeKutta(dt);
      dE = fabs(E - oldE);
  }
    
  QP[5] = Jphi/(QP[0]*QP[0]);
  while (QP[2]< 0.) QP[2] += TPi;
  while (QP[2]>TPi) QP[2] -= TPi;
}

//

OrbitIntegratorWithStats::OrbitIntegratorWithStats(Vector<double,6> StartPoint,
						   Potential *Phi,
						   const double Ttot) {
  setup(StartPoint,Phi,Ttot);
}

void  OrbitIntegratorWithStats::setup(Vector<double,6> StartPoint,
				      Potential *PhiIn,
				      const double Ttot) {
  Pot = PhiIn;
  Tmax = Ttot;
  XV_ini = StartPoint;
  Stepper.setup(XV_ini,Pot);
  //setup(StartPoint);
  Energy = Stepper.Energy();
  Lz = Stepper.AngularMomentum();
  GuidingRadius = Pot->RfromLc(Lz);
  setupDone = true;
}


void  OrbitIntegratorWithStats::setup(Vector<double,6> StartPoint) {
  XV_ini = StartPoint;
  Stepper.setup(XV_ini);
  
  Energy = Stepper.Energy();
  Lz = Stepper.AngularMomentum();
  GuidingRadius = Pot->RfromLc(Lz);
  setupDone =true;
}


// The idea is that whatever output I'm asking for, the Runge-Kutta
// integration happens here. The type tells the function what the return will
// be
int OrbitIntegratorWithStats::runGeneric(const string type,
					 Vector <double,6> *output,
					 double *tout, int N) {
  int nout_run=0;
  double t=0.;
  double outputDelt, dt=1.e-2, t_tol=1.e-4;
  Vector <double,6> XVold=XV_ini,XV;
  double tnext=0.;

  if(!setupDone) setup(XV_ini);
  
  if(Energy>0.) {
    setupDone=false;
    return -1; // Unbound - failure.
  }
  
  if(type=="NoOutput") tnext=Tmax;
  else if(type == "OutputNoTimes" || type == "OutputWithTimes") {
    
  }

  
  // initial values
  MinR = XV_ini[0];
  MaxR = XV_ini[0];
  Maxz = fabs(XV_ini[1]);
  MeanR = 0.;

  Minr = XV_ini[0]*XV_ini[0]+XV_ini[1]*XV_ini[1];
  Maxr = Minr;

  
  while(t<Tmax) {
    if(Tmax-t<dt) {
      dt = Tmax-t;
      Stepper.set_maxstep(dt);
    }
    Stepper.stepRK_by(dt);
    XV = Stepper.XV();
    double r2 = XV[0]*XV[0]+XV[1]*XV[1];
    // Could do better by more intelligent average
    MeanR += XV[0]*dt;
    // could do better with interpolation
    if(XV[0]<MinR) MinR = XV[0]; 
    if(XV[0]>MaxR) MaxR = XV[0];
    if(fabs(XV[1])> Maxz) Maxz = fabs(XV[1]);
    if(r2<Minr) Minr = r2;
    if(r2>Maxr) Maxr = r2;

    t += dt;
  }

  Maxr = sqrt(Maxr);
  Minr = sqrt(Minr);
  MeanR /= t;
  PseudoEccentricity = (Maxr-Minr)/(Maxr+Minr);

  setupDone=false;
  
  return 0; // success
}


//end of OrbitIntegrator.cc ////////////////////////////////////////////////////
