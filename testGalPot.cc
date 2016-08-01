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
    int    iso=1;
    Frequencies KNO;
    double R=8.2,z=0.014,r,dr,dz,G,ro,rho0,P,P1,P2,dPdR,dPdR2,dPdz,dPdz2,q;
    if(argc !=2) {
      cerr << "This test needs an input potential, e.g. pot/PJM16_best.Tpot\n";
      exit(1);
    }
    
    cout<< "Reading input potential described in file " << argv[1] << '\n'; 
    file.open(argv[1]);
    GalaxyPotential Phi(file);
    file.close();
 
    cout<<"Object Phi, of type `GalaxyPotential' constructed. Now testing.\n";
    cout<<"\nN.B. Everything here is in code units: kpc, Myr, M_sun\n";
    
    cout<<"\nWe work at an (approximate) solar position, R=8.2, z=0.014\n";
    cout<<" GalaxyPotential can return just the potential at that point:\n";
    cout<<" P = Phi(R,z) = ";
    P=Phi(R,z);
    cout<< P <<"\n\n";
    
    cout<< "Or the potential and derivatives (not the forces), dP/dR & dP/dz:\n";
    cout<<" P = Phi(R,z,dPdR,dPdz) = ";
    P=Phi(R,z,dPdR,dPdz);
    cout<< P <<"\n dPdR = " << dPdR << "\n dPdz = " << dPdz << "\n\n";
    
    cout<<"\nOther possible queries are density:\n rho = Phi.Density(R,z) = ";
    cout<< Phi.Density(R,z)  << "\n\n";
    
    cout<<"Circular velocity squared:\n vc2 = Phi.vcsquare(R) = ";
    cout<< Phi.vcsquare(R) << "\n\n";

    cout<<"Mass interior to radius:\n M = Phi.Mass(R) = ";
    cout<< Phi.Mass(R) << "\n\n";

    cout<<"Oort\'s constants:\n Phi.OortConstants(R,A,B)\n";
    double A,B;
    Phi.OortConstants(R,A,B);
    cout<<" A = " << A << ";\t\tB = " << B << "\n\n";

    cout<<"The radius of a circular orbit with Ang Mom Lz:\n R = Phi.RfromLc(Lz) = ";
    cout << Phi.RfromLc(2.) << " (for Lz=2.)\n\n";

    cout<<"And the inverse:\n Lz = Phi.LfromRc(R) = " << Phi.LfromRc(R) <<"\n\n";

    //    cout<<"The Laplacian of Phi, which tends towards 4 Pi G rho:\n Lap = Phi.Laplace(R,z) = "
    //	<< Phi.Laplace(R,z) << "\n\n";
    //cout << 4*Pi*Units::G*Phi.Density(R,z)<< "\n\n";
    
    cout<< "And the epicyclic frequencies (radial, vertical & azimuthal):\n";
    cout<<" KNO = Phi.KapNuOm(R) = "
	<< Phi.KapNuOm(R) << '\n'; 
}
