#include <fstream>
#include <iostream>
#include "GalPot.h"
#include <cmath>

using std::cout;
using std::cin;

extern "C" {
  GalaxyPotential* GalPot_new(char fname[]){
    ifstream file;
    file.open(fname);
    GalaxyPotential *Phi = new GalaxyPotential(file);
    file.close();
    return Phi;
  }

  void GalPot_delete(GalaxyPotential* Phi){
    delete Phi;
    return;
  }
  double GalPot_Potential_single(GalaxyPotential* Phi, double R, double z){
    double P;
    P = (*Phi)(R,z);
    return P;
  }
  void GalPot_Potential(GalaxyPotential *Phi, double P[],
  			  double R[], double z[], int length){
    for (int i=0; i<length; i++){
      P[i] = (*Phi)(R[i],z[i]);
    }
    return;
  }

  double GalPot_Potential_and_derivatives_single(GalaxyPotential* Phi,
					  double R, double z,
					  double *dPdR, double *dPdz){
    double P = (*Phi)(R,z,*dPdR,*dPdz);
    //cout << *dPdR << "\n";
    return P;
  }
  void GalPot_Potential_and_derivatives(GalaxyPotential* Phi, double P[],
					double R[], double z[],
					double dPdR[], double dPdz[],
					int length){
    for (int i=0; i<length; i++){
      P[i] = (*Phi)(R[i],z[i],dPdR[i],dPdz[i]);
    }
    return;
  }

  double GalPot_Density_single(GalaxyPotential *Phi, double R, double z){
    return (*Phi).Density(R,z);
  }
  void GalPot_Density(GalaxyPotential *Phi, double rho[], double R[], double z[],
			int length){
    for (int i=0; i<length; i++){
      rho[i] =  (*Phi).Density(R[i],z[i]);
    }
    return;
  }


  double GalPot_vcsquare_single(GalaxyPotential *Phi, double R){
    return (*Phi).vcsquare(R);
  }
  void GalPot_vcsquare(GalaxyPotential *Phi, double vc2[], double R[], int length){
    for (int i=0; i<length; i++){
      vc2[i] =  (*Phi).vcsquare(R[i]);
    }
    return;
  }

  double GalPot_Mass_single(GalaxyPotential *Phi, double R){
    return (*Phi).Mass(R);
  }
  void GalPot_Mass(GalaxyPotential *Phi, double mass[], double R[], int length){
    for (int i=0; i<length; i++){
      mass[i] =  (*Phi).Mass(R[i]);
    }
    return;
  }

  double GalPot_RfromLc_single(GalaxyPotential *Phi, double Lc){
    return (*Phi).RfromLc(Lc);
  }
  double GalPot_LfromRc_single(GalaxyPotential *Phi, double Rc){
    return (*Phi).LfromRc(Rc);
  }
  int GalPot_KapNuOm_single(GalaxyPotential *Phi, double R,
                            double*kappa, double *nu, double *Omega){
    Frequencies freqs = (*Phi).KapNuOm(R);
    *kappa = freqs(0);
    *nu = freqs(1);
    *Omega = freqs(2);
    return 0;
  }
}


int main(int argc, char const *argv[]) {
  /* code */
  char name[] = "pot/PJM17_best.Tpot";
  GalaxyPotential *Phi = GalPot_new(name);

  return 0;
}
