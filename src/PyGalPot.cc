#include <fstream>
#include <iostream>
#include "GalPot.h"
#include "OrbitIntegrator.h"
#include <cmath>

using std::cout;
using std::cin;

extern "C" {
  GalaxyPotential* GalPot_new(char fname[]){
    ifstream file;
    GalaxyPotential *Phi = NULL;
    file.open(fname);
    if(file.is_open())
      {
	Phi = new GalaxyPotential(file);
      }
    else
      {
	cerr << "Potential file "<< fname << " doesn't exist. Exiting...\n";
	exit(1);
      }
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

//-----------------------------------------------------------------------------
//
//  Routines for using OrbitIntegrator
//
//-----------------------------------------------------------------------------

  // set up a new OrbitIntegrator
  OrbitIntegratorWithStats* OrbitIntegrator_new(GalaxyPotential* Phi, double tEnd){
    Vector <double,6> XV=1.;
    OrbitIntegratorWithStats *OI = new OrbitIntegratorWithStats(XV, Phi, tEnd);
    return OI;
  }
  // delete OrbitIntegrator, freeing memory
  void OrbitIntegrator_delete(OrbitIntegratorWithStats* OI){
    delete OI;
    return;
  }

  // run the orbit integrator purely to get orbit statistics
  // Note that this can be run with multiple different stating points
  int runOrbitIntegratorforstats(OrbitIntegratorWithStats* OI, 
                                 double XV[],      // Input: R,z,phi,vR,vz,vphi
                                 double outputs[], // Output: Energy, Lz, peri, apo, Rg, e
                                 int Ndim) {       // Input: Number of input XVs
    
    Vector <double,6> XVvec;
    for (int n = 0; n<Ndim; n++) {
      for(int i=0;i<6;i++) {
        XVvec[i] = XV[6*n+i]; 
      }
      OI->setup(XVvec);
      int failure = OI->run();
      outputs[7*n+0] = OI->Energy;
      outputs[7*n+1] = OI->Lz;
      if (failure==0){
        outputs[7*n+2] = OI->Minr;
        outputs[7*n+3] = OI->Maxr;
        outputs[7*n+4] = OI->Maxz;
        outputs[7*n+5] = OI->GuidingRadius;
        outputs[7*n+6] = OI->PseudoEccentricity;
      } else {
        for(int i=2; i<7; i++) {
          outputs[7*n+i] = 0.;
        }
      }
  }
  return 0;
  }

  // run the orbit integrator purely to get orbit path(s) and statistics
  // Note that this can be run with multiple different stating points
  int runOrbitIntegratorwithPath(OrbitIntegratorWithStats* OI, 
                                  double XV[],      // Input: R,z,phi,vR,vz,vphi
                                  double tseries[], // Input: times for output
                                  int nt,           // Input: number of times
                                  double path[],    // Output: xv at t
                                  double outputs[], // Output: Energy, Lz, peri, apo, Rg, e
                                  int Ndim) {       // Input: Number of input XVs
    // I really do need that to be a 6-vector
    Vector <double,6> XVvec;
    Vector <double,6>* XV_out = new Vector <double,6>[nt];
    for (int n = 0; n<Ndim; n++) {
      for(int i=0;i<6;i++) {
        XVvec[i] = XV[6*n+i]; 
      }
      OI->setup(XVvec);
      int failure = OI->runWithOutputAtGivenTimes(XV_out,tseries,nt);
      for(int i=0; i<nt; i++) {
        for(int j=0;j<6;j++){
          if(failure == 0)
            path[n*nt*6+i*6+j] = XV_out[i][j];
          else 
            path[n*nt*6+i*6+j] = 0.;
        }
      }
      outputs[7*n+0] = OI->Energy;
      outputs[7*n+1] = OI->Lz;
      if (failure==0){
        outputs[7*n+2] = OI->Minr;
        outputs[7*n+3] = OI->Maxr;
        outputs[7*n+4] = OI->Maxz;
        outputs[7*n+5] = OI->GuidingRadius;
        outputs[7*n+6] = OI->PseudoEccentricity;
      } else {
        for(int i=2; i<7; i++) {
          outputs[7*n+i] = 0.;
        }
      }
  }
  return 0;
  }










}


int main(int argc, char const *argv[]) {
  /* code */
  char name[] = "pot/PJM17_best.Tpot";
  GalaxyPotential *Phi = GalPot_new(name);
  GalPot_delete(Phi);
  return 0;
}
