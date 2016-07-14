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
 
 if(argc<7) {
   cerr << "Input: R z phi v_R v_z v_phi\n";
   cerr << "   (distance in kpc, angle in degrees, velocity in km/s\n";
   return 1;
 }

 // Read position & velocity from input
 Vector <double,6> XV;
 for(int i=0;i!=6;i++) XV[i] = atof(argv[i+1]);

 // Convert input to code coordinates
 XV[0] *= Units::kpc;
 XV[1] *= Units::kpc;
 XV[2] *= Units::degree;
 XV[3] *= Units::kms;
 XV[4] *= Units::kms;
 XV[5] *= Units::kms;
 
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
   cout << "Input unbound in potential. Energy "
	<< OI.Energy/(Units::kms*Units::kms) << " km^2/s^2\n"
	<< std::flush;
 }


 return 0;

}
