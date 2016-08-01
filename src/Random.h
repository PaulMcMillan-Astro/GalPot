/***************************************************************************//**
\file Random.h
\brief Contains base class RandomNumberGenerator and classes Random3 Sobol
  RandomDeviate Uniform	Gaussian Exponential ExpDisk.
Basically all you need to generate random numbers.

*                                                                              *
*  Random.h                                                                    *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
*                                                                              *
*******************************************************************************/
/*                                                                             *
* class RandomNumberGenerator   base class for random number generators        *
* class Random3                 pseudo-random number generator                 *
* class Sobol                   quasi-random number generator                  *
*                                                                              *
* class RandomDeviate           base class for a random distribution           *
* class Uniform                 P(x) = 1/|b-a| for x in (a,b)
* class Gaussian                P(x) = Exp[-x^2/(2 sigma^2)]; x in [-oo,oo]    *
* class Exponential             P(x) = Exp[-x]; x in [0,oo)                    *
*                                                                              *
*******************************************************************************/

#ifndef _Random_
#define _Random_ 1

#include <algorithm>
////////////////////////////////////////////////////////////////////////////////
/** \brief    base class for random number generators	       */
class RandomNumberGenerator {
public:
    virtual double RandomDouble()        = 0;
    void   RandomDouble(double& x) { x = RandomDouble(); }
    double operator()  () 	   { return RandomDouble(); }
    void   operator()  (double& x) { RandomDouble(x); }
    //virtual ~RandomNumberGenerator();
}; 


////////////////////////////////////////////////////////////////////////////////
/** \brief base class for random distributions */
class RandomDeviate {
public:
    virtual double operator() () = 0;
    virtual double value      (const double) const = 0;
    //virtual ~RandomDeviate();
}; 

////////////////////////////////////////////////////////////////////////////////
// here are two implementation of random number generators

/** \brief pseudo-random number generator   */
class Random3 : public RandomNumberGenerator {
// pseudo-random number generator
private:
    int   *inext, *inextp;
    long  *ma;
public:
    Random3(const long);
    virtual ~Random3();
    virtual double RandomDouble();
};


/** \brief quasi-random number generator   */
class Sobol : public RandomNumberGenerator {
// quasi-random number generator
private:
    int           actl, in;
    unsigned long bits, ix, *v;
    double        fac;
    void error(const char*) const;
    void warning(const char*) const;
public:
    Sobol(const int=-1, const int=0);
    virtual ~Sobol();
    virtual double RandomDouble();
	    void   Reset       ()   { ix = in = 0; }
            int    actual      ()   { return actl; }
};

////////////////////////////////////////////////////////////////////////////////
// here are four implimentation of random distributions

/**
   \brief Random number x with P(x) = 1/|b-a| for x in (a,b). 
   Default a=0,b=1. Requires an input pointer to RandomNumberGenerator.
 */
class Uniform : public RandomDeviate {
// gives x uniformly in [a,b]
private:
    RandomNumberGenerator &r;
    double 	a,b,ba;
public:
    Uniform(RandomNumberGenerator* R, 		    // random number generator
	    const double A=0., const double B=1.)   // a,b
      : r(*R), a(std::min(A,B)), b(std::max(A,B)), ba(b-a) {}
    double lower_bound() { return a; }
    double upper_bound() { return b; }
    double operator() () { return a+ba*r(); }
    double value(const double x) const { return (a<=x && x<=b)? 1. : 0.; }
};

/**
\brief  Random number x with P(x) = Exp[-x^2/(2 sigma^2)]; x in [-oo,oo]. 
Default sigma=1.  Requires two input pointers to RandomNumberGenerator. 

*/
class Gaussian : public RandomDeviate {
// gives x in [-oo,oo] with probability proportional to exp[-x^2/(2 sigma^2)]
private:
    int    	iset;
    double 	sig, norm, gset;
    RandomNumberGenerator *R1, *R2;
public:
    Gaussian(RandomNumberGenerator*, 		// 1st random number generator
    	     RandomNumberGenerator*,  		// 2nd random number generator
	     const double=1.);			// sigma
    double operator() ();
    double sigma   () { return sig; }
    double value(const double) const;
};


/**
   \brief	Random number x with P(x) = Exp[-x/a]; x in [0,oo) 
   Default a=1. Requires an input pointer to RandomNumberGenerator.
*/
class Exponential : public RandomDeviate {
// gives x in [0,oo] with probability proportional to Exp[-x/alpha]
private:
    double alf;
    RandomNumberGenerator *Rn;
public:
    Exponential(RandomNumberGenerator* R,	// random number generator
                const double a=1.)		// alpha
	: alf(a), Rn(R) {}
    double operator() ();
    double value(const double) const;
};


#endif
