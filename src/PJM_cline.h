/***************************************************************************//**
\file PJM_cline.h
\brief Various tools for dealing with command line arguments

*                                                                              *
* PJM_cline.h                                                                  *
*                                                                              *
* C++ code written by Paul McMillan, 2011                                      *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*       Classes which output the value of a distribution function              *
*                                                                              *
*******************************************************************************/

#ifndef _PJMcline_
#define _PJMcline_ 1

#include "PJM_utils.h"


inline bool parse_comm_line(char* cl_argv,const string test,bool &out) {
  string teststr=test;
  if(string(cl_argv).substr(0,teststr.length()) == teststr) {
    char c = cl_argv[teststr.length()];
    if(c=='1' || c=='t' || c=='T' || c=='y' || c=='Y') out=true;
    else out=false;
    return true;
  }
  return false;
}

inline bool parse_comm_line(char* cl_argv,const string test,ofstream &out) {
  string teststr=test;
  if(string(cl_argv).substr(0,teststr.length()) == teststr) {
    my_open(out,string(cl_argv).substr(teststr.length()).c_str());
    return true;
  }
  return false;
}

inline bool parse_comm_line(char* cl_argv,const string test,ifstream &out) {
  string teststr=test;
  if(string(cl_argv).substr(0,teststr.length()) == teststr) {
    my_open(out,string(cl_argv).substr(teststr.length()).c_str());
    return true;
  }
  return false;
}

inline bool parse_comm_line(char* cl_argv,const string test,int &out) {
  string teststr=test;
  if(string(cl_argv).substr(0,teststr.length()) == teststr) {
    out= atoi(string(cl_argv).substr(teststr.length()).c_str());
    return true;
  }
  return false;
}

inline bool parse_comm_line(char* cl_argv,const string test,float &out) {
  string teststr=test;
  if(string(cl_argv).substr(0,teststr.length()) == teststr) {
    out= atof(string(cl_argv).substr(teststr.length()).c_str());
    return true;
  }
  return false;
}

inline bool parse_comm_line(char* cl_argv,const string test,double &out) {
  string teststr=test;
  if(string(cl_argv).substr(0,teststr.length()) == teststr) {
    out= atof(string(cl_argv).substr(teststr.length()).c_str());
    return true;
  }
  return false;
}

inline bool parse_comm_line(char* cl_argv,const string test,string &out) {
  string teststr=test;
  if(string(cl_argv).substr(0,teststr.length()) == teststr) {
    out=string(cl_argv).substr(teststr.length());
    return true;
  }
  return false;
}



inline string my_get_rootname(char* input) {
  string tmpstr=string(input); 
  int tmpi;
  
  for(int i=tmpstr.size();tmpstr[i]!='.' && i>0;i--) tmpi=i-1;
  
  if(tmpi) return tmpstr.substr(0,tmpi);
  else     return tmpstr;
}

inline int split_comma_val_input(char* input) {
  string str = string(input);
  size_t comma1 = str.find(','), comma2;
  int nout=0;
  cerr << str.substr(0,comma1) << ",";
  nout++;
  while(comma1 != string::npos) {
    comma2 = str.find(',',comma1+1);  
    cerr << str.substr(comma1+1,comma2-comma1-1).c_str() << ","; 
    comma1 = comma2;
    nout++;
  }
  cerr << "\n";
  return nout;
}

inline int split_comma_val_input(char* input, double* table, const int ntabmax)
{
  if(ntabmax<1) return 0;
  string str = string(input);
  size_t comma1 = str.find(','), comma2;
  int nout=0;
  table[nout++] = atof(str.substr(0,comma1).c_str());

  while(comma1 != string::npos && nout != ntabmax) {
    comma2 = str.find(',',comma1+1);  
    table[nout++] = atof(str.substr(comma1+1,comma2-comma1-1).c_str()); 
    comma1 = comma2;
  }
  if(comma1 != string::npos) { 
    cerr << "More comma separated command line arguments than capacity\n";
  }
  return nout;
}

inline int split_comma_val_input(char* input, float* table, const int ntabmax)
{
  if(ntabmax<1) return 0;
  string str = string(input);
  size_t comma1 = str.find(','), comma2;
  int nout=0;
  table[nout++] = atof(str.substr(0,comma1).c_str());

  while(comma1 != string::npos && nout != ntabmax) {
    comma2 = str.find(',',comma1+1);  
    table[nout++] = atof(str.substr(comma1+1,comma2-comma1-1).c_str()); 
    comma1 = comma2;
  }
  if(comma1 != string::npos) { 
    cerr << "More comma separated command line arguments than capacity\n";
  }
  return nout;
}

inline int split_comma_val_input(char* input, int* table, const int ntabmax)
{
  if(ntabmax<1) return 0;
  string str = string(input);
  size_t comma1 = str.find(','), comma2;
  int nout=0;
  table[nout++] = atoi(str.substr(0,comma1).c_str());

  while(comma1 != string::npos && nout != ntabmax) {
    comma2 = str.find(',',comma1+1);  
    table[nout++] = atoi(str.substr(comma1+1,comma2-comma1-1).c_str()); 
    comma1 = comma2;
  }
  if(comma1 != string::npos) { 
    cerr << "More comma separated command line arguments than capacity\n";
  }
  return nout;
}

#endif
