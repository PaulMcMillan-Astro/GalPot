/*******************************************************************************
*                                                                              *
*  findManyOrbitPropertiesfromGalactic.cc                                      *
*                                                                              *
*  C++ code written by Paul McMillan, 2007-                                    *
*  Lund Observatory, Lund University.                                          *
*  address: Box 43, SE-221 00 Lund, Sweden                                     *
*  e-mail:  paul@astro.lu.se                                                   *
*                                                                              *
*******************************************************************************/


#include "GalPot.h"
#include "PJMCoords.h"
#include "OrbitIntegrator.h"

#include <ctime>

int main(int argc,char *argv[])  
{

  ifstream file,data;
  ofstream output;
  string potfile = "pot/PJM16_best.Tpot",line;
  Potential *Phi;

  OmniCoords OC;
  

  file.open(potfile.c_str());
  if(!file){
    cerr << "Input file does not exist. Filename: " << potfile << "\n";
    return 0;
  }
  Phi = new GalaxyPotential(file);
  file.close();
  
  if(argc<3) {
    cerr << "Input: input_file output_file\n";
    cerr << "Input is an ascii file giving position and motion in galactic coordinates:\n"
	 << "\tdistance l b v_los mu_l* mu_b\n";
    cerr << "\t(distance in kpc, angles in degrees, velocity in km/s, proper motion in mas/yr)\n";
    cerr << " output is in kpc (distances), km^2/s^2 (energy), kpc km/s (angular momentum)\n";
    return 0;
  }

  data.open(argv[1]);
  if(!file){
    cerr << "Input file does not exist. Filename: " << string(argv[1]) << "\n";
    return 0;
  }

  output.open(argv[2]);

  output << "#MinR MaxR Maxz Minr Maxr MeanR Energy AngMom\n" << std::flush;
  
  Vector <double,6> XV=5.;

  Vector <double,6> GalacticCoords;
  
  OrbitIntegratorWithStats OI(XV, Phi, 10000.);

  
  
  while(getline(data,line)) {
    if(line[0] != '#') {
      std::stringstream ss(line);
      for(int i=0;i!=6;i++) ss >> GalacticCoords[i];

      GalacticCoords[0] *= Units::kpc;
      GalacticCoords[1] *= Units::degree;
      GalacticCoords[2] *= Units::degree;
      GalacticCoords[3] *= Units::kms;
      GalacticCoords[4] *= Units::masyr;
      GalacticCoords[5] *= Units::masyr;

      XV = OC.GCYfromHGP(GalacticCoords);

      // One could instead use:
      // OC.takeHGP_units(GalacticCoords); // (without converting units)
      // XV = OC.giveGCY();
      
      OI.setup(XV);
      if(OI.run() == 0) {
	output << OI.MinR << ' ' << OI.MaxR<< ' ' << OI.Maxz<< ' '
	       << OI.Minr << ' ' << OI.Maxr<< ' ' << OI.MeanR<< ' '
	       << OI.Energy/(Units::kms*Units::kms)<< ' '
	       << OI.Lz/(Units::kms*Units::kpc) << '\n';
      }
      
    }

  }




}
