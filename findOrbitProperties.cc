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


#include "falPot.h"
#include "OrbitIntegrator.h"

#include <ctime>

int main(int argc,char *argv[])  
{

 ifstream file;
 string potfile = "pot/PJM16_best.Tpot";
 Potential *Phi;
 

 file.open(potfile.c_str());
 if(!file) {
   cerr << "Input file does not exist. ";
   cerr << "Filename: " << potfile << "\n";
   return 0;
 }
 Phi = new GalaxyPotential(file);
 file.close();
 
 if(argc<7) {
   cerr << "Input: R z phi v_R v_z v_phi\n";
   cerr << "   (distance in kpc, angle in degrees, velocity in km/s\n";
   return 0;
 }
 
 Vector <double,6> XV;
 for(int i=0;i!=6;i++) XV[i] = atof(argv[i+1]);
 
 XV[0] *= Units::kpc;
 XV[1] *= Units::kpc;
 XV[2] *= Units::degree;
 XV[3] *= Units::kms;
 XV[4] *= Units::kms;
 XV[5] *= Units::kms;
 

 OrbitIntegratorWithStats OI(XV, Phi, 10000.);
 clock_t start = clock();
 int Ntot = 1000;

 //for(int i=0;i!=Ntot;i++) { OI.run(); }
 
 //cerr << (double)(clock()-start)/((double)(Ntot*CLOCKS_PER_SEC)) <<'\n';
			    
   
 if(OI.run() == 0) {
   cerr << "Guiding Centre radius: " << OI.GuidingRadius << '\n';
   cerr << "Minimum, Maximum Cylindrical radius: "<< OI.MinR << ','
	<< OI.MaxR << '\n'
	<< "Maximum z: " << OI.Maxz << '\n'
	<< "Minimum, Maximum Spherical radius: "<< OI.Minr << ','
	<< OI.Maxr << '\n'
	<< "Energy: " << OI.Energy/(Units::kms*Units::kms) << " km^2/s^2\n"
	<< "Angular Momentum (about symmetry axis): "
	<< OI.Lz/(Units::kms*Units::kpc)
	<< " kpc km/s\n"
       	<< "Mean Cylindrical radius: "<< OI.MeanR << '\n';
 }




}
