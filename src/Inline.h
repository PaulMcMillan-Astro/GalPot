//-----------------------------------------------------------------------------+
//                                                                             |
// Inline.h                                                                    |
//                                                                             |
// C++ code                                                                    |
//                                                                             |
// written by Walter Dehnen, 1994-2001                                         |
// e-mail:   dehnen@mpia-hd.mpg.de                                             |
// address:  Max-Planck Institut für Astronomie, Königstuhl 69117 Heidelberg   |
//           Germany                                                           |
//                                                                             |
//-----------------------------------------------------------------------------+
//
// Contains lots of stuff I'm not going to look at too closely. Hopefully 
// nothing too dangerous...
//

#ifndef _Inline_def_
#define _Inline_def_ 1

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

using std::cerr;
////////////////////////////////////////////////////////////////////////////////
#ifndef ebug
void Numerics_error(const char*) __attribute__ ((noreturn));
#else
void Numerics_error(const char*);
#endif

inline void Numerics_error(const char* msgs)
{
    cerr << " ERROR in Numerics: " << msgs << '\n';
#ifndef ebug
    exit(1);
#endif
}

////////////////////////////////////////////////////////////////////////////////
inline int Numerics_message(const char* msgs)
{
  cerr << " WARNING in Numerics: " << msgs << '\n';
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>

//using namespace std;

using std::ofstream;
using std::ifstream;
using std::fstream;
using std::istream;
using std::ostream;
using std::ios;
using std::string;

//------------------------------------------------------------------------------
#define TM template<class S> inline

//------------------------------------------------------------------------------
TM S    abs     (const S x)              { return (x<0)? -x : x; }
TM int  sign    (const S x)              { return (x<0)? -1:((x>0)? 1:0 ); }
TM S    sign    (const S x, const S s)   { return ( s>0 )?  abs(x) : -abs(x); }
//TM S    min     (const S x, const S y)   { return (x<y)? x : y; }
//TM S    max     (const S x, const S y)   { return (x>y)? x : y; }
TM S    mod     (const S x, const S y)   { return x-y*int(x/y); }
TM S    square  (const S x)              { return x*x; }
TM S    cube    (const S x)              { return x*x*x; }
TM void swap    (S&a, S&b)               { register S t=a; a=b; b=t; }
TM S    pow_uint(const S x, const unsigned int n)
{
  if(n==0) return S(1);
  register S z=x, y=(n%2)? x : S(1);
  for(register unsigned int i=n>>1; i; i>>=1) { z*=z; if(i%2) y*=z; }
  return y;
}
TM S    pow_int (const S x, const int i)
{
  if(i>=0) return pow_uint(x,i);
#ifdef C_AND_F
  if(x==0) { printf("pow_int(): negative power of zero"); exit(1); }
#else // C_AND_F
  if(x==0) { cerr<<"pow_int(): negative power of zero"; exit(1); }
#endif // C_AND_F
  return pow_uint(S(1)/x,-i);
}

#undef TM


//------------------------------------------------------------------------------
#ifndef C_AND_F
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void SwallowRestofLine(istream& from)
{
    char c;
    do from.get(c); while( from.good() && c !='\n');
}
//------------------------------------------------------------------------------
inline const char* stndrdth(const int i)
{
    register int ia= (i<0)? -i:i;
    switch( ia % 100 ) {
        case 11: 
        case 12: 
        case 13: return "th";
        default:
        switch( ia % 10 ) {
            case 1:  return "st";
            case 2:  return "nd"; 
            case 3:  return "rd";
            default: return "th";
        }   
    }
}
//------------------------------------------------------------------------------
#endif // C_AND_F


template<class S> inline S sqr(const S x)              { return square(x); }
template<class S> inline S pow(const S x, const int i) { return pow_int(x,i); }



////////////////////////////////////////////////////////////////////////////////
#endif
