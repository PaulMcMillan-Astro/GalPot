/*******************************************************************************
*                                                                              *
*  findManyOrbitPropertiesfromEquatorialwErrors.cc                             *
*                                                                              *
*  C++ code written by Paul McMillan, 2007-                                    *
*  Lund Observatory, Lund University.                                          *
*  address: Box 43, SE-221 00 Lund, Sweden                                     *
*  e-mail:  paul@astro.lu.se                                                   *
*                                                                              *
*******************************************************************************/


#include "GalPot.h"
#include "Random.h"
#include "PJMCoords.h"
#include "OrbitIntegrator.h"

#include <ctime>

#include <algorithm>    // std::sort
#include <vector>       // std::vector

void OutputMedianAndUpperLower(ofstream &output, std::vector<double> &VecIn) {
  std::sort(VecIn.begin(),VecIn.end());
  int n_tests = VecIn.size();
  double median;
  if(n_tests%2 == 1) median =  VecIn[n_tests/2];
  else               median =  0.5*(VecIn[n_tests/2]+VecIn[n_tests/2]);
  output << median << ' ';
  //15.87, 84.13 percentiles
  if(n_tests<4) output << VecIn[n_tests-1]-median << ' ' << median-VecIn[0] << ' ';
  
  for(int i=0;i!=2;i++) {
    double percentile = (i==0)? 0.1587 : 0.8413;
    double place = percentile*n_tests - 0.5;
    int iplace = int(place);
    output << fabs(VecIn[iplace] + (place-iplace)*(VecIn[iplace+1]-VecIn[iplace]) - median) << ' ';
   
  }
}

int main(int argc,char *argv[])  
{

  ifstream file,data;
  ofstream output;
  int n_tests;
  string potfile = "pot/PJM16_best.Tpot",line;
  Potential *Phi;

  // Start some random number generators
  // see Random.h, Random.cc
  time_t tm = time(NULL);
  Random3 R3(int(tm)+3), R3b(int(tm)+678);
  Gaussian GaussianRandom(&R3,&R3b);

    
  OmniCoords OC;
  

  file.open(potfile.c_str());
  if(!file){
    cerr << "Input file does not exist. Filename: " << potfile << "\n";
    return 0;
  }
  Phi = new GalaxyPotential(file);
  file.close();
  
  if(argc<4) {
    cerr << "Input: input_file n_tests output_file\n";
    cerr << "input_file is an ascii file giving position and motion in equatorial coordinates, with uncertainties. n_tests is number of Monte Carlo sample used per star\n"
	 << "\tdistance distance _err RA RA_err DEC DEC_err v_los v_los_err mu_a* mu_a*_err mu_d mu_d_err\n";
    cerr << "\t(distance in kpc, angles in degrees, velocity in km/s, proper motion in mas/yr)\n";
    cerr << " Output is in kpc (distances), km^2/s^2 (energy), kpc km/s (angular momentum)\n";
    cerr << " Quoted values are medians from Monte Carlo sample, with minus and plus 1sigma values found from 15.87 and 84.13 percentiles.\n";
    return 0;
  }

  data.open(argv[1]);
  if(!file){
    cerr << "Input file does not exist. Filename: " << string(argv[1]) << "\n";
    return 0;
  }

  n_tests = atoi(argv[2]);
  output.open(argv[3]);

  output << "#MinR MinRMinus MinRPlus MaxR MaxRMinus MaxRPlus Maxz MaxzMinus MaxzPlus Minr MinrMinus MinrPlus Maxr MaxrMinus MaxrPlus MeanR MeanRMinus MeanRPlus Energy EnergyMinus EnergyPlus AngMom AngMomMinus AngMomPlus \n" << std::flush;

  // Tables to hold results
  std::vector<double> MinR, MaxR, Maxz, Minr, Maxr, MeanR, Energy, AngMom;
  
  Vector <double,6> XV=1.;

  Vector <double,6> EquatorialCoords, EquatorialCoordsErr, EquatorialCoordsTmp;
  
  OrbitIntegratorWithStats OI(XV, Phi, 10000.);

  
  while(getline(data,line)) {
    int DangerPoint = n_tests-(int(0.8413*n_tests - 0.5)+1);
    int BadPoints=0;
    if(line[0] != '#') {
      std::stringstream ss(line);
      for(int i=0;i!=6;i++) ss >> EquatorialCoords[i] >> EquatorialCoordsErr[i];

      for(int i=0;i!=n_tests;i++) {
	// Add uncertainties
	for(int j=0;j!=6;j++)
	  EquatorialCoordsTmp[j] = EquatorialCoords[j]
	    + EquatorialCoordsErr[j] * GaussianRandom();
	
	while(EquatorialCoordsTmp[0] < 0.) {
	  // I'm not letting you have a negative distance. Don't be silly.
	   EquatorialCoordsTmp[0] = EquatorialCoords[0]
	     + EquatorialCoordsErr[0] * GaussianRandom();
	}
	EquatorialCoordsTmp[0] *= Units::kpc;
	EquatorialCoordsTmp[1] *= Units::degree;
	EquatorialCoordsTmp[2] *= Units::degree;
	EquatorialCoordsTmp[3] *= Units::kms;
	EquatorialCoordsTmp[4] *= Units::masyr;
	EquatorialCoordsTmp[5] *= Units::masyr;

	XV = OC.GCYfromHEQ(EquatorialCoordsTmp);
	
	OI.setup(XV);
	if(OI.run() == 0) {
	  MinR.push_back(OI.MinR/Units::kpc); 
	  MaxR.push_back(OI.MaxR/Units::kpc);
	  Maxz.push_back(OI.Maxz/Units::kpc);
	  Minr.push_back(OI.Minr/Units::kpc);
	  Maxr.push_back(OI.Maxr/Units::kpc);
	  MeanR.push_back(OI.MeanR/Units::kpc);
	  Energy.push_back(OI.Energy/(Units::kms*Units::kms));
	  AngMom.push_back(OI.Lz/(Units::kpc*Units::kms)); 
	} else {
	  BadPoints++;
	  // Orbit unbound. Put in some sensible values
	  MinR.push_back(XV[0]/Units::kpc); 
	  MaxR.push_back(1.e10);
	  Maxz.push_back(1.e10);
	  Minr.push_back(sqrt(XV[0]+XV[0]+XV[1]*XV[1])/Units::kpc);
	  Maxr.push_back(1.e10);
	  MeanR.push_back(1.e10);
	  Energy.push_back(OI.Energy/(Units::kms*Units::kms));
	  AngMom.push_back(OI.Lz/(Units::kpc*Units::kms)); 
	}
      }
      // when finished
      // Output median and pm 1sigma equivalent percentiles

      if(BadPoints>=DangerPoint) {
	std::cerr << "WARNING: Output unreliable\t"
		  << "Too many unbound orbits for star with coordinates "
		  << EquatorialCoords << "\n";
      }
      
      OutputMedianAndUpperLower(output,MinR);
      OutputMedianAndUpperLower(output,MaxR);
      OutputMedianAndUpperLower(output,Maxz);
      OutputMedianAndUpperLower(output,Minr);
      OutputMedianAndUpperLower(output,Maxr);
      OutputMedianAndUpperLower(output,MeanR);
      OutputMedianAndUpperLower(output,Energy);
      OutputMedianAndUpperLower(output,AngMom);
      output << '\n';
      
      
    }
  }



}
