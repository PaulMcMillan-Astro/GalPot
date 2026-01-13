/*******************************************************************************
*                                                                              *
*  findManyOrbitProperties.cc                                                  *
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

int main(int argc,char *argv[])  
{

  ifstream file,data;
  ofstream output;
  string potfile = "pot/PJM16_best.Tpot",line;
  Potential *Phi;
 

 // Read potential from file 
  file.open(potfile.c_str());
  if(!file){
    cerr << "Input file does not exist. Filename: " << potfile << "\n";
    return 1;
  }
  Phi = new GalaxyPotential(file);
  file.close();
  
  if(argc<3) {
    cerr << "Input: input_file output_file\n";
    cerr << "input_file must contain columns: R z v_R v_z v_phi\n";
    cerr << "   (distance in kpc, velocity in km/s)\n";
    cerr << " output is in kpc (distances), km^2/s^2 (energy), kpc km/s (angular momentum)\n";
    cerr << " Warning: If energy is negative, other properties will be incorrect";
    return 1;
  }

  // Open file for input
  data.open(argv[1]);
  if(!file){
    cerr << "Input file does not exist. Filename: " << string(argv[1]) << "\n";
    return 0;
  }

  // Open for output
  output.open(argv[2]);

  // Write header 
  output << "#MinR MaxR Maxz Minr Maxr MeanR Energy AngMom\n" << std::flush;

  // Setup class
  Vector <double,6> XV=1.;
  OrbitIntegratorWithStats OI(XV, Phi, 10000.);

  int nline;
  // read file
  while(getline(data,line)) {
    // Skip commented lines (which start with #)
    if(line[0] != '#') {
      std::stringstream ss(line);
      
      for(int i=0;i!=5;i++)
	if(i<2) ss >> XV[i];
	else ss >> XV[i+1];
       // Convert input to code coordinates
      XV[0] *= Units::kpc;
      XV[1] *= Units::kpc;
      //XV[2] *= Units::degree;
      XV[3] *= Units::kms;
      XV[4] *= Units::kms;
      XV[5] *= Units::kms;
      
      OI.setup(XV);
      // run integration
      if(OI.run() == 0) {
      // Output results if successful
	output << OI.MinR<< ' ' << OI.MaxR<< ' ' << OI.Maxz<< ' '
	       << OI.Minr<< ' ' << OI.Maxr<< ' ' << OI.MeanR<< ' '
	       << OI.Energy/(Units::kms*Units::kms)<< ' '
	       << OI.Lz/(Units::kms*Units::kpc) << '\n' << std::flush;
      } else {
	cerr << "Failure for line: " << line
	     << "\nEnergy=" << OI.Energy/(Units::kms*Units::kms) << '\n';
      }
    }

  }


  return 0;

}
