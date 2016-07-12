/*******************************************************************************
*                                                                              *
*  testfalPot.cc                                                               *
*                                                                              *
*  C++ code written by Walter Dehnen, 1995-96,                                 *
*                      Paul McMillan, 2007-08,                                 *
*  Oxford University, Department of Physics, Theoretical Physics.              *
*  address: 1 Keble Road, Oxford OX1 34P, United Kingdom                       *
*  e-mail:  p.mcmillan1@physics.ox.ac.uk                                       *
*                                                                              *
*******************************************************************************/

#include <fstream>
#include <iostream>
#include "falPot.h"
#include <cmath>

using std::cout;
using std::cin;


int main(int argc,char *argv[])  
{
    ifstream file;
    int    iso=1;
    Frequencies KNO;
    double R=8.,z=0.007,r,dr,dz,G,ro,rho0,P,P1,P2,dPdR,dPdR2,dPdz,dPdz2,q;
    if(argc !=2) {
      cerr << "This test needs an input potential, e.g. pot/DB97Mod1.Tpot\n";
      exit(1);
    }
    
    cout<< "Reading input potential described in file " << argv[1] << '\n'; 
    file.open(argv[1]);
    GalaxyPotential Phi(file);
    file.close();
 
    cout<<" Object Phi, of type `GalaxyPotential' constructed. Now testing.\n";
    cout<<"\n We work at an (approximate) solar position, R=8.0, z=0.007\n";
    cout<<" GalaxyPotential can return just the potential at that point\n";
    cout<<" P = Phi(R,z) = ";
    P=Phi(R,z);
    cout<< P <<"\n";
    cout<< "Or the potential and derivatives (not the forces), dP/dR & dP/dz\n";
    cout<<" P = Phi(R,z,dPdR,dPdz) = ";
    P=Phi(R,z,dPdR,dPdz);
    cout<< P <<"\n dPdR = " << dPdR << "\n dPdz = " << dPdz << '\n';
    cout<<"\n Other possible queries are\n rho = Phi.Density(R,z)\n";
    cout<<" vc2 = vcsquare(R)\n M=Phi.Mass(R)\n ";
    cout<<" Phi.OortConstants(R,A,B), where A & B";
    cout<<" have their usual meanings in this context\n";
    cout<<" R = Phi.RfromLc(Lz), ";
    cout<<"the radius of a circular orbit with Ang Mom Lz\n";
    cout<<" Lz = Phi.LfromRc(R), the inverse\n";
    cout<<" Lap = Phi.Laplace(R,z), gives Laplace(Phi), tends towards 4 Pi G\n";
    cout<<" KNO = Phi.KapNuOm(R), epicycle frequencies, radial, vertical & azimuthal\n"; 
}
