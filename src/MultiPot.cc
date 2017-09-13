// Combining potentials.

#include "MultiPot.h"


double MultiPotential::operator() (const double R, const double z) const
{
  double Phi=0;

  for(int i=0;i!=npot;i++)
    Phi += (*(PotList[i]))(R,z);

  return Phi;

}


double MultiPotential::operator() (const double R, const double z, double &dPdR, double &dPdz) const
{
  double Phi=0, dPdR_tmp, dPdz_tmp;
  dPdR = 0.;
  dPdz = 0.;

  for(int i=0;i!=npot;i++) {
    Phi += (*(PotList[i]))(R,z, dPdR_tmp,dPdz_tmp);
    dPdR += dPdR_tmp;
    dPdz += dPdz_tmp;
  }

  return Phi;
}



// double MultiPotential::RfromLc(const double L_in, double* dR) const
// {
//   bool more=false;
//   double R,lR=0.,dlR=0.001,z,dPR,dPz,P,LcR,oldL,L=fabs(L_in);
//   R=exp(lR);
//   P= (*this)(R,0.,dPR,dPz);
//   LcR=sqrt(R*R*R*dPR);
//   if(LcR == L) return R;
//   if(L>LcR) more=true;
//   oldL=LcR;

//   for( ; ; ) {
//     lR += (more)? dlR : -dlR;
//     R=exp(lR);
//     P= (*this)(R,0.,dPR,dPz);
//     LcR=sqrt(R*R*R*dPR);
//     if(LcR == L) return R;
//     if((L< LcR && L>oldL) ||(L>LcR && L<oldL)){
// 	R=(more)? exp(lR-0.5*dlR) : exp(lR+0.5*dlR);
// 	return R;}
//     oldL=LcR;
//   }

// }

Frequencies MultiPotential::KapNuOm(const double R) const {
  Frequencies KNO = 0., KNOtmp;
  for(int i=0;i!=npot;i++) {
    KNOtmp = (*(PotList[i])).KapNuOm(R);
    for(int j=0;j!=3;j++) KNO[j] += KNOtmp[j]*KNOtmp[j];
  }

  for(int j=0;j!=3;j++) KNO[j] = sqrt(KNO[j]);

  return KNO;

}
