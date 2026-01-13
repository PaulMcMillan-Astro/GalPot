/***************************************************************************//**
*                                                                              *
* MiyamotoNagaiPot.cc                                                          *
*                                                                              *
* C++ code written by Paul McMillan, 2006-                                     *
*                                                                              *
*******************************************************************************/
#include "MiyamotoNagaiPot.h"
#include <cmath>

double MiyamotoNagaiPotential::operator() (const double R, const double z) const
{
  double AZB=A+sqrt(z*z+Bq);
  return -GM/sqrt(R*R+AZB*AZB);
}

double MiyamotoNagaiPotential::operator() (const double R, const double z, double& dPdR,
			       double& dPdz) const
{
  double  ZB=sqrt(z*z+Bq), AZB = A + ZB,
    F = 1./(R*R+AZB*AZB),rtF=sqrt(F);
  dPdR = GM*R*F*rtF;
  dPdz = ZB? GM*z*F*rtF*AZB/ZB : 0.;
  return -GM*rtF;
}


Frequencies MiyamotoNagaiPotential::KapNuOm(const double R) const
{
  double dPR,dPz,P;
  P = (*this)(R,0.,dPR,dPz);
  double F = 1./(R*R+ABq),rtF=sqrt(F),
    om2    = dPR/R, 
    nu2    = GM*ABoB*F*rtF,
    kappa2 = GM*(F*rtF - 3.*R*R*F*F*rtF) + 3.*om2;
  Frequencies output= sqrt(kappa2);
  output[1] = sqrt(nu2);
  output[2] = sqrt(om2);

  return output;
}

ostream& operator<< (ostream& to, const MiyamotoNagaiPotential& P)
{
  to << "Miyamoto-Nagai potential with mass=" << P.GM/Units::G
     << " a=" << P.A << " b=" << sqrt(P.Bq) << '\n';
    return to;
}








