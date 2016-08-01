/*******************************************************************************
*                                                                              *
*  Random.cc                                                                   *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  dehnen@thphys.ox.ac.uk                                              *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <cmath>
#include "Random.h"
#include "Numerics.h"


////////////////////////////////////////////////////////////////////////////////
// class Random3 *************************************************************//
////////////////////////////////////////////////////////////////////////////////

const long   mbig  = 1000000000,
	     mseed = 161803398;
const double fac   = 1./double(mbig);

Random3::Random3(const long idum)
{
    inext  = new int;
    inextp = new int;
    ma     = new long[56];

    register long  mj,mk;
    register int   i,ii,k;

    mj     = mseed - (idum<0 ? -idum : idum);
    mj    %= mbig;
    ma[55] = mj;
    mk     = 1;
    for(i=1; i<=54; i++) {
	ii     = (21*i) % 55;
	ma[ii] = mk;
	mk     = mj-mk;
	mj     = ma[ii];
	if(mk<0) mk += mbig;
    }
    for(k=1; k<=4; k++)
	for(i=1; i<=55; i++) {
	    ma[i] -= ma[1+(i+30) % 55];
	    if(ma[i]<0) ma[i] += mbig;
	}
    *inext  = 0;
    *inextp = 31;
}

Random3::~Random3()
{
    delete inext;
    delete inextp;
    delete[] ma;
}

double Random3::RandomDouble ()
{
    register long   mj;
    register double r;
    do {
	if(++(*inext)  >= 56) *inext  = 1;
	if(++(*inextp) >= 56) *inextp = 1;
	mj = ma[*inext] - ma[*inextp];
	while(mj<0) mj+=mbig;
	ma[*inext] = mj;
	r=mj*fac;
    } while (r<0. || r>1.);
    return r;
}

////////////////////////////////////////////////////////////////////////////////
// class Sobol ***************************************************************//
////////////////////////////////////////////////////////////////////////////////

static int  setb  = 30;
const  int  MO    = 52;
static char f[MO] ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
const  int  d[MO] ={1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,
		    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
		    8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
const  int  p[MO] ={0,1,1,2,1,4,2,4,7,11,13,14,1,13,16,19,22,25,
		    1,4,7,8,14,19,21,28,31,32,37,41,42,50,55,56,59,62,
		    14,21,22,38,47,49,50,52,56,67,70,84,97,103,115,122};

//
// with this installation not more than 52 objects (Sobol sequences) can be
// created simultaneously.
//

void Sobol::error(const char* msgs) const
{
    cerr << " error in Sobol: " << msgs << '\n';
#ifndef ebug
    exit(1);
#endif
}

void Sobol::warning(const char* msgs) const
{
    cerr << " warning in Sobol: " << msgs << '\n';
}

void set_bits(const int BITS)
{
    setb = (BITS<=0) ? 30 : BITS;
}

Sobol::Sobol(const int ACTL, const int BITS)
{
// set degree and polynomial from the static tables. It's important, that
// no other object with the same degree and polynomial is currently existing,
// because its quasi random numbers would be the same as ours. If ACTL on 
// input is set >=0 we don't care and construct the sequence No. ACTL.
    if(ACTL>=0 && ACTL < MO)
	actl=ACTL;
    else {
        for(actl=0; actl<MO && f[actl]; actl++) {}
        if(actl>=MO)  error("trying to create the 53th object");
    }
    f[actl] += 1;
    if(BITS==0)
	bits = setb;
    else if((bits=BITS)<10)
	warning("creating object with less than 10 bits");
    in       = 0;
    ix       = 0;
    fac      = 1./(1L<<bits);
    int degs = d[actl];
    int poly = p[actl];
// seed initial Mi; i=1,...,degs
// these must be odd integer numbers less than 2^i.
// Finally the direction numbers are Vi = 2^(bits-i) * Mi
    register int i,i2,ip,l;
    register unsigned long vi;
    v = (new unsigned long[bits])-1;
    for(i=1,i2=2; i<=degs; i++,i2<<=1) {
	if(i2<=poly) 
	    vi = 1;
	else {
	    vi = i2 - poly;
	    if(!(vi&1)) vi-= 1;
	}
	v[i] = vi << (bits-i);
    }
// now use the recurrence (Press et al. 1992, eq. 7.7.2) to create
// the remaining direction numbers. With Vi = 2^(bits-i) Mi it reads
// V[i] = (a[1]*V[i-1]) XOR (a[2]*V[i-2]) XOR ... XOR (a[q-1]*V[i-q+1])
//	  XOR ( V[i-q] XOR V[i-q]/2^q )
    for(i=degs+1; i<=int(bits); i++) {
	ip = poly;
	vi = v[i-degs];
	vi^= (vi>>degs);
	for(l=degs-1; l>0; l--) {
	    if(ip&1) vi ^= v[i-l];
	    ip >>= 1;
	}
	v[i] = vi;
    }
}

Sobol::~Sobol()
{
    delete[] (v+1);
    f[actl]  = 0;
}

double Sobol::RandomDouble ()
{
    register unsigned long im=in++, j;
    for(j=1; j<=bits; j++) {
	if( !(im&1) ) break;
	im >>= 1;
    }
    if(j>bits) error("trying to call more than 2^BITS times");
    ix^= v[j];
    return double(ix)*fac;
}

////////////////////////////////////////////////////////////////////////////////
// class Gaussian ************************************************************//
////////////////////////////////////////////////////////////////////////////////

Gaussian::Gaussian(RandomNumberGenerator* r1,    // 1st random number generator
		   RandomNumberGenerator* r2,    // 2nd random number generator
		   const double s)               // sigma
: sig(s), R1(r1), R2(r2)
{ 
    iset = 0;  
    norm = 1./(sqrt(TPi)*sig);
}

double Gaussian::operator() ()
{
    if(iset) {
	iset = 0;
	return gset;
    } else {
	register double v1,v2,rsq,fac;
	do {
	    v1  = 2 * R1->RandomDouble() - 1.;
	    v2  = 2 * R2->RandomDouble() - 1.;
	    rsq = v1*v1 + v2*v2;
	} while (rsq>=1. || rsq <=0. );
	fac  = sig*sqrt(-2.*log(rsq)/rsq);
	gset = v1*fac;
	iset = 1;
	return v2*fac;
    }
}

double Gaussian::value(const double x) const 
{
    return norm * exp(-0.5*pow(x/sig,2));
}

////////////////////////////////////////////////////////////////////////////////
// class Exponential *********************************************************//
////////////////////////////////////////////////////////////////////////////////

double Exponential::operator() ()
{
    return -alf * log( Rn->RandomDouble() );
}

double Exponential::value(const double x) const
{
    return exp(-x/alf)/alf;
}
