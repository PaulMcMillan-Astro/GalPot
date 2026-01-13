/*******************************************************************************
*                                                                              *
*  findManyOrbitPropertiesfromEquatorial.cc                                    *
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
    cerr << "Input is an ascii file giving position and motion in equatorial coordinates (J2000):\n"
	 << "\tdistance RA DEC v_los mu_a* mu_d\n";
    cerr << "\t(distance in kpc, angles in degrees, velocity in km/s, proper motion in mas/yr)\n";
    cerr << " output is in kpc (distances), km^2/s^2 (energy), kpc km/s (angular momentum)\n";
    //cerr << "   (distance in kpc, angle in degrees, velocity in km/s\n";
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

  Vector <double,6> EquatorialCoords;
  
  OrbitIntegratorWithStats OI(XV, Phi, 10000.);

  
  
  while(getline(data,line)) {
    if(line[0] != '#') {
      std::stringstream ss(line);
      for(int i=0;i!=6;i++) ss >> EquatorialCoords[i];
      
      EquatorialCoords[0] *= Units::kpc;
      EquatorialCoords[1] *= Units::degree;
      EquatorialCoords[2] *= Units::degree;
      EquatorialCoords[3] *= Units::kms;
      EquatorialCoords[4] *= Units::masyr;
      EquatorialCoords[5] *= Units::masyr;

      XV = OC.GCYfromHEQ(EquatorialCoords);
      
      OI.setup(XV);
      if(OI.run() == 0) {
	output << OI.MinR<< ' ' << OI.MaxR<< ' ' << OI.Maxz<< ' '
	       << OI.Minr<< ' ' << OI.Maxr<< ' ' << OI.MeanR<< ' '
	       << OI.Energy/(Units::kms*Units::kms)<< ' '
	       << OI.Lz/(Units::kms*Units::kpc) << '\n';
      }
      
    }

  }




}
