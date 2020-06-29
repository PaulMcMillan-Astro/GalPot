/*******************************************************************************
*                                                                              *
* C++ code written by Paul McMillan, 2014-                                     *
*                                                                              *
*******************************************************************************/
#include "KeplerPot.h"
#include <cmath>

double KeplerPotential::operator() (const double R, const double z) const
{
  double rsq=R*R+z*z;
  return -GM/(sqrt(rsq));
}

double KeplerPotential::operator() (const double R, const double z, double& dPdR,
			       double& dPdz) const
{
  double rsq = R*R+z*z,
    r = sqrt(rsq), ir = 1./r,
    dPdr = GM/rsq;
  dPdR = dPdr*R/r;
  dPdz = dPdr*z/r;
  return -GM*ir;

}


Frequencies KeplerPotential::KapNuOm(const double R) const
{
  
  double tmp = sqrt(GM/(R*R*R));

  Frequencies epi = tmp; // all three components equal
  return epi;

}
