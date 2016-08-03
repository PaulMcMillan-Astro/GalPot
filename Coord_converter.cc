/***************************************************************************//**
\file Coord_converter.cc
\brief Converts between coordinates (various Heliocentric & Galactocentric systems)

*                                                                              *
*  Coord_converter.cc                                                          *
*                                                                              *
*  C++ code written by Paul McMillan 2011-,                                    *
*                      Walter Dehnen, 1995-96,                                 *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/GalPot                        *
*                                                                              *
*******************************************************************************/
//
// v 0.1 Coord_converter created    15/9/11
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "PJM_cline.h"
using std::cout;


#include "PJMCoords.h"

int main(int argc,char *argv[])
{

  if(argc<4) { 
    cerr << "Inputs: XXXfromXXX input_file output_file (nskip=0) "
	 << "(code_units_input=f) (code_units_output=f) " 
	 << "(R0="<< GalactoConstants::Rsun <<") "
	 << "(z0="<< GalactoConstants::zsun <<") "
	 << "(v0="<< GalactoConstants::vcsun*Units::kms_i<<")"
	 << " (epoch=2000)\n" 
	 << argv[0] << " h for help\n"; 
    if(argc>1) if(argv[1][0] == 'h') {
      cerr << "Program to convert between different coordinate systems (see "
	   << "PJMCoords.h)\n"
	   << "XXX refer to coordinate systems. Possible choices:\n"
	   << "\tGCA\t\tGalactocentric CArtesian (x,y,z,vx,vy,vz)\n"
	   << "\tGCY\t\tGalactocentric CYlindrical (R,z,phi,vR,vz,vphi)\n"
	   << "\tLSR\t\tLocal Standard of Rest (x,y,z,vx,vy,vz)\n"
	   << "\tHCA\t\tHeliocentric CArtesian (x,y,z,u,v,w)\n"
	   << "\tHGP\t\tHeliocentric Galactic Polar (r,l,b,vr,mu_l*,mu_b)\n"
	   << "\tHEQ\t\tHeliocentric EQuatorial polar (r,a,d,vr,mu_a*,mu_d)\n\n"
	   << "Input file must contain coordinates in each row, with nskip "
	   << "values before the first coord value in each case\n"
	   << "Output is just the new coordinates\n"
	   << "code_units_(in/out)put true if (in/out)put values "
	   << "in kpc, Myr, radians\n"
	   << "R0,z0 given in kpc, v0 in km/s (and negative),\n"
	   << "Non-code units are kpc for distances, degrees for angles,"
	   << "km/s for velocities, mas/yr for proper motions\n";
    }
    exit(0);
  }
  bool code_in=false, code_out=false; 
  int n_in=1, infile_offset=0, n_per_line=0, n_excess=0, codeunits=0;
  double R0 = GalactoConstants::Rsun, 
    v0 = GalactoConstants::vcsun*Units::kms_i,
    z0 = GalactoConstants::zsun,
    original_epoch, epoch=0.; // not known here

  for(int i=4;i<argc;i++) {
    bool understood=false;
    if(parse_comm_line(argv[i],"nskip=",infile_offset))        understood=true;
    if(parse_comm_line(argv[i],"code_units_input=",code_in))   understood=true;
    if(parse_comm_line(argv[i],"code_units_output=",code_out)) understood=true;
    if(parse_comm_line(argv[i],"R0=",R0))                      understood=true;
    if(parse_comm_line(argv[i],"v0=",v0))                      understood=true;
    if(parse_comm_line(argv[i],"z0=",z0))                      understood=true;
    if(parse_comm_line(argv[i],"epoch=",epoch))                understood=true;
    if(!understood) {
      cerr << "Input "<<argv[i]<<" not understood\n"; 
      exit(0);
    }
  }

  double line_pre[1000], line_post[1000];
  OmniCoords OC;
  OC.give_epoch(original_epoch);
  OC.change_sol_pos(R0,z0);   
  OC.change_vc(v0*Units::kms);
  if(epoch!=0. && epoch != original_epoch) OC.change_epoch(epoch);
  // OC.change_vsol(10*Units::kms, 10*Units::kms, 10*Units::kms);
  // left here in case you need to change the peculiar motion of the Sun.

  
  vec6 input, output;
  ifstream from;
  ofstream to;
  my_open(from,argv[2]);
  my_open(to,argv[3]);
  n_in = how_many_lines(from);
  n_per_line = entrys_in_line(from);
  n_excess = n_per_line-infile_offset-6;
  if(n_excess<0) {
      cerr << "too few columns\n"; exit(1);
  }
  string str = string(argv[1]);
  if(str.length() != 10) {
    cerr << "requirement "<< argv[1] << "not understood";
    exit(1);
  }
  string intype = str.substr(7,3);
  string outtype = str.substr(0,3);

  for(int i=0;i!=n_in;i++) {
    for(int j=0;j<infile_offset;j++) from >> line_pre[j];
    for(int j=0;j!=6;j++) from >> input[j];
    for(int j=0;j!=n_excess;j++) from >> line_post[j];
  
    if(code_in) {
      if(intype=="HEQ") OC.take_HEQ(input);
      else if(intype=="HGP") OC.take_HGP(input);
      else if(intype=="HCA") OC.take_HCA(input);
      else if(intype=="LSR") OC.take_LSR(input);
      else if(intype=="GCA") OC.take_GCA(input);
      else if(intype=="GCY") OC.take_GCY(input);
      else {
	cerr << "requirement "<< argv[1] << "not understood\n";
	exit(1);
      }
    } else {
      if(intype=="HEQ") OC.take_HEQ_units(input);
      else if(intype=="HGP") OC.take_HGP_units(input);
      else if(intype=="HCA") OC.take_HCA_units(input);
      else if(intype=="LSR") OC.take_LSR_units(input);
      else if(intype=="GCA") OC.take_GCA_units(input);
      else if(intype=="GCY") OC.take_GCY_units(input);
      else {
	cerr << "requirement "<< argv[1] << "not understood\n";
	exit(1);
      }
    }
    //for(int j=0;j<infile_offset;j++) to << line_pre[j] << ' ';
    if(!code_out) {
      if(outtype=="HEQ") to << OC.give_HEQ_units()<< " ";
      else if(outtype=="HGP") to << OC.give_HGP_units()<< " ";
      else if(outtype=="HCA") to << OC.give_HCA_units()<< " ";
      else if(outtype=="LSR") to << OC.give_LSR_units()<< " ";
      else if(outtype=="GCA") to << OC.give_GCA_units()<< " ";
      else if(outtype=="GCY") to << OC.give_GCY_units()<< " ";
      else {
	cerr << "requirement "<< argv[1] << "not understood\n";
	exit(1);
      }
    } else {
    
      if(outtype=="HEQ") to << OC.give_HEQ()<< " ";
      else if(outtype=="HGP") to << OC.give_HGP()<< " ";
      else if(outtype=="HCA") to << OC.give_HCA()<< " ";
      else if(outtype=="LSR") to << OC.give_LSR()<< " ";
      else if(outtype=="GCA") to << OC.give_GCA()<< " ";
      else if(outtype=="GCY") to << OC.give_GCY()<< " ";
      else {
	cerr << "requirement "<< argv[1] << "not understood\n";
	exit(1);
      }
    }
    //for(int j=0;j!=n_excess;j++) to << line_post[j] << ' ';
    to << '\n' << std::flush;
  }

}
