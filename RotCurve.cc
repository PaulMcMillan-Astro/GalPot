/*******************************************************************************
*                                                                              *
*  testGalPot.cc                                                               *
*                                                                              *
*  C++ code written by Walter Dehnen, 1995-96,                                 *
*                      Paul McMillan, 2007-,                                   *
*  Lund Observatory, Lund University.                                          *
*  address: Box 43, SE-221 00 Lund, Sweden                                     *
*  e-mail:  paul@astro.lu.se                                                   *
*                                                                              *
*******************************************************************************/

#include <fstream>
#include <iostream>
#include "GalPot.h"

using std::cout;
using std::cin;


int main(int argc,char *argv[])
{
    ifstream file;
    ofstream outfile;

    int    npoints=100;
    double Rmin=0.1, Rmax=10;
    
    if(argc < 3) {
      cerr << "\nInput: potfile outfile (Rmin Rmax npoints)\n";
      cerr << "e.g. pot/PJM16_best.pot out.tab 0.1 10 100\n\n";
      std::cerr << "Last three input are optional."
                <<" If missing, default values (0.1 10 100) used" << "\n\n";
      exit(1);
    }

    cout<< "Reading input potential described in file " << argv[1] << '\n';
    file.open(argv[1]);
    GalaxyPotential Phi(file);
    file.close();
    std::cout << "Writing to " << argv[2] << '\n';
    outfile.open(argv[2]);
    if(argc>3) Rmin = atof(argv[3]);
    if(argc>4) Rmax = atof(argv[4]);
    if(argc>5) npoints = atoi(argv[5]);

    std::cout << "For "<< npoints << " points between "
              << Rmin << " and " << Rmax << '\n';

    double dR = (npoints>1)? (Rmax-Rmin)/(npoints-1) : (Rmax-Rmin);


    for (int i = 0; i < npoints; i++) {
    /* code */
      double R = Rmin + i*dR;
      outfile << R << ' ' << sqrt(Phi.vcsquare(R))/Units::kms << "\n";
    }

  std::cout << "Done." << '\n';
}
