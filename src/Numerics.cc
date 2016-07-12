/*******************************************************************************
*                                                                              *
* Numerics.cc                                                                  *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
* Lund Observatory, Lund University.                                           *
* address: 1 Keble Road, Oxford, OX1 3NP, United Kingdom.                      *
* e-mail:  w.dehnen1@physics.ox.ac.uk                                          *
*                                                                              *
*******************************************************************************/

#include <cmath>
#include "Numerics.h"
#include "FreeMemory.h"



////////////////////////////////////////////////////////////////////////////////
typedef int*     Pint;
typedef float*   Pflt;
typedef float**  PPflt;
typedef double*  Pdbl;
typedef double** PPdbl;


////////////////////////////////////////////////////////////////////////////////
double qbulir(double(*func)(double), const double a, const double b, 
	      const double eps_, double& err)
/*------------------------------------------------------------------------------
Quadrature program using the Bulirsch sequence and rational extrapolation. The
algorithm is puplished in Bulirsch & Stoer, Num. Math. 9, 271-278 (1967), where
a routine in ALGOL is given. This routine is a straightforward translation into
C++.
CAUTION: 
Do not use this routine for integrating low order polynomials (up to fourth
order) or periodic functions with period equal to the interval of integration
or linear combinations of both.
INPUT:  func   pointer to function to be integrated.
        a,b    lower and upper boundaries of the integration interval;
        eps    desired relativ accuracy;
OUTPUT: return approximated value for the integral;
        err    actual relative error of the return value.
------------------------------------------------------------------------------*/
{
    register double ba=b-a;
    if(ba==0.) return 0.;

    register int    i,n=2,nn=3,mx=25,m,mr, bo,bu=0,odd=1;
    register double c,d1,ddt,den,e,eps,eta=1.e-7,gr,hm,nt,
		    sm,t,t1,t2,t2a,ta,tab=0.,tb,v=0.,w;
    double          d[7],dt[7];

    while(eta+1. != 1.) eta *=0.5;
    eta  *=2.;                       // eta = actual computing accuracy

    eps   = max(eps_,eta);
    sm    = 0.;
    gr    = 0.;
    t1    = 0.;
    t2    = 0.5*((*func)(a)+(*func)(b));
    t2a   = t2;
    tb    = abs(t2a);
    c     = t2*ba;
    dt[0] = c;

    for(m=1;m<=mx;m++) {           // iterate over the refinements
	bo = (m>=7);
	hm = ba/n;
	if(odd) {
	    for(i=1;i<=n;i+=2) {
		w  = (*func)(a+i*hm);
		t2+= w;
		tb+= abs(w);
	    }
	    nt  = t2;
	    tab = tb * abs(hm);
	    d[1]=16./9.;
	    d[3]=64./9.;
	    d[5]=256./9.;
	} else {
	    for(i=1;i<=n;i+=6) {
		w  = i*hm;
		t1+= (*func)(a+w) + (*func)(b-w);
	    }
	    nt  = t1+t2a;
	    t2a =t2;
	    d[1]=9./4.;
	    d[3]=9.;
	    d[5]=36.;
	}
	ddt  =dt[0];
	t    =nt*hm;
	dt[0]=t;
	nt   =dt[0];
	if(bo) {
	    mr  =6;
	    d[6]=64.;
	    w   =144.;
	} else {
	    mr  =m;
	    d[m]=n*n;
	    w   =d[m];
	}
	for(i=1;i<=mr;i++) {
	    d1 =d[i]*ddt;
	    den=d1-nt;
	    e  =nt-ddt;
	    if(den != 0.) {
		e /= den;
		v  = nt*e;
		nt = d1*e;
		t += v;
	    } else {
		nt = 0.;
		v  = 0.;
	    }
	    ddt  = dt[i];
	    dt[i]= v;
	}
	ta = c;
	c  = t;
	if(!bo) t -= v;
	v  = t-ta;
	t += v;
	err= abs(v);
	if(ta<t) {
	    d1 = ta;
	    ta = t;
	    t  = d1;
        }
	bo = bo || (ta<gr && t>sm);
	if(bu && bo && err < eps*tab*w) break;
	gr = ta;
	sm = t;
	odd= !odd;
	i  = n;
	n  = nn;
	nn = i+i;
	bu = bo;
	d[2]=4.;
        d[4]=16.;
    }
    v = tab*eta;
    if(err<v) err = v;
    if(m==mx) Numerics_error("qbulir exceeding maximum of iterations");
    return c;
}
////////////////////////////////////////////////////////////////////////////////
void GaussLegendre(Pdbl x, Pdbl w, const int n)
{
    register double eps;
    for(eps=1.e-10; (eps+1.)!=1.; eps*=0.5);
    eps  =1.e-10;                       // eps != actual computing accuracy
                                  // because that was crashing the f'ing thing
    register int j,i,m=(n+1)/2;
    register double z1,z,pp,p3,p2,p1;
    for (i=0;i<m;i++) {
	z=cos(Pi*(i+0.75)/(n+0.5));
	do {
	    p1 = 1.0;
	    p2 = 0.0;
	    for(j=0;j<n;j++) {
		p3 = p2;
		p2 = p1;
		p1 = ( (2*j+1)*z*p2 - j*p3 ) / double(j+1);
	    }
	    pp = n * (z*p1-p2) / (z*z-1.0);
	    z1 = z;
	    z  = z1 - p1 / pp;
	} while (abs(z-z1)>eps);
	x[i]     =-z;
	x[n-1-i] = z;
	w[i]     = 2. / ((1.0-z*z)*pp*pp);
	w[n-1-i] = w[i];
    }
}
////////////////////////////////////////////////////////////////////////////////
void LegendrePeven(double* p, const double x, const int np)
// based on a routine from J.J. Binney
// evaluates even Legendre Polys up to l=2*(np-1) at x
{
    register int    n,l,l2;
    register double x2=x*x;
    p[0] = 1.;
    p[1] = 1.5*x2-0.5;
    for(n=2; n<np; n++) {
	l = 2*(n-1);
	l2= 2*l;
	p[n] = - p[n-2] * l*(l-1)         / double((l2+1)*(l2-1))
	       + p[n-1] * (x2-(l2*l+l2-1) / double((l2-1)*(l2+3)));
        p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
}
////////////////////////////////////////////////////////////////////////////////
void dLegendrePeven(double* p, double* d, const double x, const int np)
// based on a routine from J.J. Binney
// evaluates even Legendre Polys and its derivs up to l=2*(np-1) at x
{
    register int    n,l,l2;
    register double x2=x*x;
    p[0] = 1.;
    d[0] = 0.;
    p[1] = 1.5*x2-0.5;
    d[1] = 1.5;
    for(n=2; n<np; n++) {
	l = 2*(n-1);
	l2= 2*l;
	p[n] = - p[n-2] * l*(l-1)         / double((l2+1)*(l2-1))
	       + p[n-1] * (x2-(l2*l+l2-1) / double((l2-1)*(l2+3)));
        p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
        d[n] = - d[n-2] * l*(l-1)         / double((l2+1)*(l2-1))
               + d[n-1] * (x2-(l2*l+l2-1) / double((l2-1)*(l2+3)))
	       + p[n-1];
        d[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
    x2 = 2*x;
    for(n=0; n<np; n++)
	d[n] *= x2;
}


////end of Numerics.cc//////////////////////////////////////////////////////////
