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

void OutputMedianAndUpperLower(ofstream &output, std::vector<double> &VecIn, bool trailingcomma=true) {
  std::sort(VecIn.begin(),VecIn.end());
  int n_tests = VecIn.size();
  double median;
  if(n_tests%2 == 1) median =  VecIn[n_tests/2];
  else               median =  0.5*(VecIn[n_tests/2]+VecIn[n_tests/2]);
  output << median << ',';
  //15.87, 84.13 percentiles
  if(n_tests<4) {
      output << VecIn[n_tests-1]-median << ',' << median-VecIn[0] ;
      if(trailingcomma) output << ',';
  }
  for(int i=0;i!=2;i++) {
    double percentile = (i==0)? 0.1587 : 0.8413;
    double place = percentile*n_tests - 0.5;
    int iplace = int(place);
    output << fabs(VecIn[iplace] + (place-iplace)*(VecIn[iplace+1]-VecIn[iplace]) - median);
    if(i==0 || trailingcomma) output << ',';

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

  //output << "#MinR MinRMinus MinRPlus MaxR MaxRMinus MaxRPlus Maxz MaxzMinus MaxzPlus Minr MinrMinus MinrPlus Maxr MaxrMinus MaxrPlus MeanR MeanRMinus MeanRPlus Energy EnergyMinus EnergyPlus AngMom AngMomMinus AngMomPlus Eccentricity EccentricityPlus EccentricityMinus\n" << std::flush;
  output << "X,XPlus,XMinus,Y,YPlus,YMinus,Z,ZPlus,ZMinus,vX,vXPlus,vXMinus,vY,vYPlus,vYMinus,vZ,vZPlus,vZMinus,"
         << "R,RPlus,RMinus,z,zPlus,zMinus,vR,vRPlus,vRMinus,vz,vzPlus,vzMinus,vphi,vphiPlus,vphiMinus\n";
  // Tables to hold results
  std::vector<double> X,Y,Z,vX,vY,vZ,R,z,vR,vz,vphi;
  //std::vector<double> MinR, MaxR, Maxz, Minr, Maxr, MeanR, Energy, AngMom;
  //std::vector<double> Ecc;

  Vector <double,6> XV=1.,XV_cart=1.;

  Vector <double,6> EquatorialCoords, EquatorialCoordsErr, EquatorialCoordsTmp;



  while(getline(data,line)) {
    int DangerPoint = n_tests-(int(0.8413*n_tests - 0.5)+1);
    int BadPoints=0;
    X.clear();   Y.clear();  Z.clear();
    vX.clear(); vY.clear(); vZ.clear();
    R.clear();  z.clear();
    vR.clear(); vz.clear(); vphi.clear();
    if(line[0] != '#') {
      std::stringstream ss(line);
      for(int i=0;i!=6;i++) ss >> EquatorialCoords[i] >> EquatorialCoordsErr[i];
      //cerr << EquatorialCoords << '\n';
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
    XV_cart = OC.HCAfromHEQ(EquatorialCoordsTmp);
	XV = OC.GCYfromHEQ(EquatorialCoordsTmp);
    X.push_back(XV_cart[0]/Units::kpc);
    Y.push_back(XV_cart[1]/Units::kpc);
    Z.push_back(XV_cart[2]/Units::kpc);
    vX.push_back(XV_cart[3]/Units::kms);
    vY.push_back(XV_cart[4]/Units::kms);
    vZ.push_back(XV_cart[5]/Units::kms);

    R.push_back(XV[0]/Units::kpc);
    z.push_back(XV[1]/Units::kpc);
    vR.push_back(XV[3]/Units::kms);
    vz.push_back(XV[4]/Units::kms);
    vphi.push_back(XV[5]/Units::kms);
    }
      // when finished
      // Output median and pm 1sigma equivalent percentiles

      if(BadPoints>=DangerPoint) {
	    std::cerr << "WARNING: Output unreliable\t"
		          << "Too many unbound orbits for star with coordinates "
		          << EquatorialCoords << "\n";
      }
      OutputMedianAndUpperLower(output,X);
      OutputMedianAndUpperLower(output,Y);
      OutputMedianAndUpperLower(output,Z);
      OutputMedianAndUpperLower(output,vX);
      OutputMedianAndUpperLower(output,vY);
      OutputMedianAndUpperLower(output,vZ);
      OutputMedianAndUpperLower(output,R);
      OutputMedianAndUpperLower(output,z);
      OutputMedianAndUpperLower(output,vR);
      OutputMedianAndUpperLower(output,vz);
      OutputMedianAndUpperLower(output,vphi,false);

      output << '\n'<< std::flush;


    }
  }



}
