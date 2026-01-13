/*******************************************************************************
*                                                                              *
*  findOrbitProperties.cc                                                      *
*                                                                              *
*  C++ code written by Paul McMillan, 2007-                                    *
*  Lund Observatory, Lund University.                                          *
*  address: Box 43, SE-221 00 Lund, Sweden                                     *
*  e-mail:  paul@astro.lu.se                                                   *
*                                                                              *
*******************************************************************************/


#include "GalPot.h"
#include "OrbitIntegrator.h"

#include <ctime>

using std::cout;

int main(int argc,char *argv[])
{

 ifstream file;
 string potfile = "pot/PJM16_best.Tpot";
 Potential *Phi;

 // Read potential from file
 file.open(potfile.c_str());
 if(!file) {
   cerr << "Input file does not exist. ";
   cerr << "Filename: " << potfile << "\n";
   return 1;
 }
 Phi = new GalaxyPotential(file);
 file.close();

 if(argc<6) {
   cerr << "Input: R z v_R v_z v_phi\n";
   cerr << "   (distance in kpc, velocity in km/s)\n";
   return 1;
 }

 // Read position & velocity from input (In galactocentric cylindrical coordinates)
 Vector <double,6> XV=0.;

 // Convert input to code coordinates
 XV[0] = atof(argv[1]) * Units::kpc; // R
 XV[1] = atof(argv[2]) * Units::kpc; // z
 XV[2] = 0.* Units::degree; // phi coordinate irrelevant
 XV[3] = atof(argv[3]) * Units::kms; // v_R
 XV[4] = atof(argv[4]) * Units::kms; // v_z
 XV[5] = atof(argv[5]) * Units::kms; // v_phi

 // Set up Integrator class
 OrbitIntegratorWithStats OI(XV, Phi, 10000.);

 // run integration
 int IntegrationFail = OI.run();

 // Output results
 if(IntegrationFail == 0) {
   cout << "Guiding Centre radius: " << OI.GuidingRadius << '\n' << std::flush;
   cout << "Minimum, Maximum Cylindrical radius: "<< OI.MinR << ','
	<< OI.MaxR << '\n'
	<< "Maximum z: " << OI.Maxz << '\n'
	<< "Minimum, Maximum Spherical radius: "<< OI.Minr << ','
	<< OI.Maxr << '\n'
	<< "Energy: " << OI.Energy/(Units::kms*Units::kms) << " km^2/s^2\n"
	<< "Angular Momentum (about symmetry axis): "
	<< OI.Lz/(Units::kms*Units::kpc)
	<< " kpc km/s\n"
       	<< "Mean Cylindrical radius: "<< OI.MeanR << '\n' << std::flush;
 } else {
	// This line should never be reached (change 2026)
   cout << "Input unbound in potential. Energy "
	<< OI.Energy/(Units::kms*Units::kms) << " km^2/s^2\n"
	<< std::flush;
 }


 return 0;

}
