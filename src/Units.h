/*******************************************************************************
*                                                                              *
*  Units.h                                                                     *
*                                                                              *
* C++ code written by Walter Dehnen, 1995/96,                                  *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                           *
*                                                                              *
********************************************************************************
*                                                                              *
* The unit system employed has:                                                *
*                                                                              *
* angle:      measured in radian,                                              *
* length:     measured in kilo parsec,                                         *
* time:       measured in Mega years,                                          *
* mass:       measured in solar masses.                                        *
*                                                                              *
* This implies the following dimensions                                        *
*                                                                              *
* quantity        dimension / seize                using other units           *
*------------------------------------------------------------------------------*
* angular vel.  1 Myr^-1	                  = 977.775320024919 km/s/kpc  *
* velocity      1 kpc/Myr                         = 977.775320024919 km/s      *
* action/mass   1 kpc^2/Myr                       = 977.775320024919 kpc*km/s  *
* potential     1 kpc^2/Myr^2                     = 956044.576449833 (km/s)^2  *
* acceleration  1 kpc/Myr^2                                                    *
* G             4.49865897 E-12 kpc^3/Myr^2/Msun                               *
* 4 Pi G        5.65318158 E-11 kpc^3/Myr^2/Msun                               *
*                                                                              *
* Note that this implies handy numbers for typical quantities in the Galaxy:   *
*                                                                              *
* velocities 		are of the order of 0.1 to 0.4                         *
* radii      		are of the order of 0.1 to 100                         *
* dynamical times	are of the order of 1   to 1000                        *
* G times total mass    is  of the order of 1                                  *
*                                                                              *
********************************************************************************
*                                                                              *
* ALL the relevant quantities are declared as public members of class Units,   *
* which is constructed without arguments. After the definition of class Units  *
* below, a const object                                                        *
*	const Units U;                                                         *
* is declared at the end of this file. This should be used by any user, e.g.   *
* U.kms refers to km/s in the basic units (1 km/s = 0.0010227299). This whole  *
* construction (i) avoids global variables (apart from U) and hence potential  *
* problems with the global name space, and (ii) ensures that nobody can change *
* the values of the quantities (since U is const). Additionally, the construc- *
* tion of further objects of type Units is disabled and will create a warning  *
* on compilation and an error at run time.                                     *
*                                                                              *
*******************************************************************************/

#ifndef _Units_set_
#define _Units_set_ 1

#include "Pi.h"
#include <string>

using std::string;

class Units {
public:
// names of the basic units
    const string angle_unit, length_unit, time_unit, mass_unit;
// basic units:
    const double radian, kpc, Myr, Msun;
// other units of the basic dimensions:
    const double rad,degree,arcmin,arcsec,mas,anghr,angmin,angsec; // angle
    const double cm,meter,km,ly,pc,Mpc,AU;			// length
    const double sec,hour,day,yr,hyr,century,Gyr;		// time
    const double gram,kg;					// mass
// derived units:
    const double kpcMyr,kpcGyr,AUyr,kms,c_light; 		// velocity
    const double radMyr,kmskpc,masyr,ashyr,secyr,asyr; 		// angle vel.
    const double pc2,kpc2,cm2;					// area
    const double pc3,kpc3,cm3;					// volume
// physical constants:
//  const double c_light; (already defined above)		// vel. of light
    const double G,Grav,fPiG;					// Newtons G
// inverse quantities:
    const double rad_i,radian_i,degree_i,arcmin_i,arcsec_i,mas_i,
	   anghr_i,angmin_i,angsec_i;
    const double cm_i,meter_i,km_i,AU_i,ly_i,pc_i,kpc_i,Mpc_i;
    const double sec_i,hour_i,day_i,yr_i,hyr_i,century_i,Myr_i,Gyr_i;
    const double gram_i,kg_i,Msun_i;
    const double kpcMyr_i,kpcGyr_i,AUyr_i,kms_i,c_light_i;
    const double radMyr_i,kmskpc_i,masyr_i,ashyr_i,secyr_i,asyr_i;
    const double pc2_i,kpc2_i,cm2_i;
    const double pc3_i,kpc3_i,cm3_i;
    const double G_i,Grav_i,fPiG_i;
	
// constructor takes no arguments
    Units() :
    // names of basic units
	angle_unit    	("radian"),
	length_unit   	("kpc"),
	time_unit     	("Myr"),
	mass_unit     	("M_sun"),
    // basic units
	radian		(1.),
	kpc		(1.),
	Myr		(1.),
	Msun		(1.),
    // angle
	rad             ( radian),			// radian
	degree		( rad * TPi / 360.),		// degrees (360=circle)
	arcmin		( degree / 60.),		// arc minutes
	arcsec		( arcmin / 60.),		// arc seconds
	mas		( 0.001 * arcsec),		// milli arc seconds
	anghr		( rad * TPi / 24.),		// angle hour(24=circle)
	angmin		( anghr / 60.),			// anlge minutes
	angsec		( angmin / 60.),		// angle seconds
    // length
	cm		( kpc * 3.240778828894144e-22),	// centimeter
	meter		( 1.e2 * cm),			// meter
	km		( 1.e5 * cm),			// kilo meter
	ly              ( 3.0660669364447e-4 * kpc),	// light year
	pc		( 1.e-3 * kpc),			// parsec
	Mpc		( 1.e3 * kpc),			// Mega parsec
	AU		( arcsec * pc),			// astronomical unit
    // time
	sec		( 3.168753556551954e-14 * Myr),	// time second
	hour		( 3600  * sec),			// time hour
	day		( 24    * sec),			// time day
	yr		( 1.e-6 * Myr),			// year
	hyr		( 1.e-4 * Myr),			// hundred years
	century 	( hyr),				// century
	Gyr		( 1.e3  * Myr),			// giga year
    // mass
	gram		( Msun / 1.989e33),		// gram
	kg		( 1.e3 * gram),			// kilogram
    // velocity
	kpcMyr		( kpc / Myr),			// kpc per Myr
	kpcGyr		( kpc / Gyr),			// kpc per Gyr
	AUyr		( AU / yr),			// AU per yr
	kms		( km / sec),			// km per second
	c_light		( ly / yr),			// speed of light
    // angle velocity
	radMyr		( radian / Myr),		// radian per Myr
	kmskpc		( kms / kpc),			// kms per kpc
	masyr		( mas / yr),			// milli arcsec per year
	ashyr		( arcsec / hyr),		// arcsec per century
	secyr		( angsec / yr),			// angsec per yr
	asyr		( arcsec / yr),			// arcsec per yr
    // area
	pc2		( pc * pc),			// square parsec
	kpc2		( kpc * kpc),			// square kilo parsec
	cm2		( cm * cm),			// square centimeter
    // volume
	pc3		( pc2 * pc),			// cubic parsec
	kpc3		( kpc2 * kpc),			// cubic kilo parsec
	cm3		( cm2 * cm),			// cubic centimeter
    // constant of gravity
	G      		( 4.498658966346282e-12),	// Newtons G
	Grav   		( G),				// Newtons G
	fPiG   		( 5.653181583871732e-11),	// 4 Pi G
    // inverse quantities (useful for transformations)
    // inverse angle
	rad_i   	( 1./rad),
	radian_i   	( 1./radian),
	degree_i 	( 1./degree),
	arcmin_i  	( 1./arcmin),
	arcsec_i  	( 1./arcsec),
	mas_i  		( 1./mas),
	anghr_i   	( 1./anghr),
	angmin_i  	( 1./angmin),
	angsec_i  	( 1./angsec),
    // inverse length
	cm_i      	( 1./cm),
	meter_i		( 1./meter),
	km_i      	( 1./km),
	AU_i      	( 1./AU),
	ly_i		( 1./ly),
	pc_i      	( 1./pc),
	kpc_i     	( 1./kpc),
	Mpc_i     	( 1./Mpc),
    // inverse time
	sec_i     	( 1./sec),
	hour_i		( 1./hour),
	day_i		( 1./day),
	yr_i      	( 1./yr),
	hyr_i     	( 1./hyr),
	century_i	( 1./century),
	Myr_i     	( 1./Myr),
	Gyr_i     	( 1./Gyr),
    // inverse mass
	gram_i    	( 1./gram),
	kg_i		( 1./kg),
	Msun_i    	( 1./Msun),
    // inverse velocity
	kpcMyr_i	( 1./kpcMyr),
	kpcGyr_i	( 1./kpcGyr),
	AUyr_i		( 1./AUyr),
	kms_i    	( 1./kms),
	c_light_i  	( 1./c_light),
    // inverse angle velocity
        radMyr_i        ( 1./radMyr),
	kmskpc_i   	( 1./kmskpc),
	masyr_i		( 1./masyr),
	ashyr_i		( 1./ashyr),
	secyr_i		( 1./secyr),
	asyr_i		( 1./asyr),
    // inverse area
	pc2_i		( 1./pc2),
	kpc2_i		( 1./kpc2),
	cm2_i		( 1./cm2),
    // inverse volume
	pc3_i		( 1./pc3),
	kpc3_i		( 1./kpc3),
	cm3_i		( 1./cm3),
    // inverse of constant of gravity
	G_i      	( 1./G),
	Grav_i   	( G_i),
	fPiG_i   	( 1./fPiG)
    {}
};

static Units U;

#endif
