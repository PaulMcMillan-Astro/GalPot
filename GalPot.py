

#  Code written by Paul McMillan 2016-
from ctypes import *
import numpy as np
import os.path
from collections import namedtuple

#import os
#print(os.path.dirname(os.path.realpath(__file__)))

lib = cdll.LoadLibrary(os.path.dirname(os.path.realpath(__file__)) +
                       '/obj/libPyGalPot.so')



class GalaxyPotential:

    """Interface for Dehnen's GalaxyPotential C++ code

    Parameters:
    ----------
    filename : string
        Name of the parameter file for the GalaxyPotential.

    This is put here as a (hopefully) convenient method of using GalPot 
    in this new Pythonic world.
    The original C++ code is still the main element.

    Everything is given in GalaxyPotential's internal units, with are 
    kpc, Myr and M_solar

    Methods
    ----------
    Potential(R,z) :            returns P
        Potential at postion(s)
    Potential_derivatives(R,z): returns P, dPdR, dPdz 
        Potential and derivatives at position(s)
    ForceRz(R,z) :              returns f_R,f_z 
        specific force from potential at position(s)
    Acceleration(xv_cyl) :      returns a_R, a_z, a_phi 
        acceleration of moving body in potential
    XVDerivative(xv_cyl):       returns dxv/dt 
        derivative of phase space position (for orbit integration)
    Density(R,z):               returns rho 
        density at position(s)
    Vcirc2(R):                  returns vcirc^2 
        square of circular velocity at R (single value or array)
    Vcirc(R):                   returns vcirc
        circular velocity at R (single value or array)
    KapNuOm(R)                  returns kappa,nu,omega 
        epicycle frequencies at R (single value only)
    RfromLc(Lz)                 returns R
        radius of circular orbit with angular momentum Lz
    LfromRc(R)                  returns Lz 
        angular momentum of circular orbit with radius R

    Members
    -----------
    kpc_Myr_to_km_s: conversion factor for velocities (from internal units to km/s)

    Notes:
    Thanks to Til Piffl, from whom I took some of this code.

    """

    kpc_Myr_to_km_s = 977.77

    def __init__(self,filename):
        """Intialise interface for Dehnen's GalaxyPotential C++ code

        Parameters:
            filename: file containing parameters of galaxy potential
        """

        # Declare various things from the C++ code
        lib.GalPot_new.restype = c_void_p
        lib.GalPot_new.argtypes = [c_char_p]

        lib.GalPot_delete.restype = c_void_p
        lib.GalPot_delete.argtypes = [c_void_p]
        lib.GalPot_Potential_single.restype = c_double
        lib.GalPot_Potential_single.argtypes = [c_void_p, c_double,c_double]

        lib.GalPot_Potential_and_derivatives_single.restype = c_double
        lib.GalPot_Potential_and_derivatives_single.argtypes = [c_void_p,c_double,c_double,
                                                                POINTER(c_double),POINTER(c_double)]

        lib.GalPot_Density_single.restype = c_double
        lib.GalPot_Density_single.argtypes = [c_void_p,c_double,c_double]

        lib.GalPot_vcsquare_single.restype = c_double
        lib.GalPot_vcsquare_single.argtypes = [c_void_p,c_double]

        lib.GalPot_Mass_single.restype = c_double
        lib.GalPot_Mass_single.argtypes = [c_void_p,c_double]

        lib.GalPot_Potential.restype = c_double
        lib.GalPot_Potential.argtypes = [c_void_p, POINTER(c_double),POINTER(c_double),POINTER(c_double),c_int]

        lib.GalPot_Potential_and_derivatives.restype = c_double
        lib.GalPot_Potential_and_derivatives.argtypes = [c_void_p,POINTER(c_double),POINTER(c_double),
                                                        POINTER(c_double),POINTER(c_double),
                                                        POINTER(c_double),c_int]

        lib.GalPot_Density.restype = c_double
        lib.GalPot_Density.argtypes = [c_void_p,POINTER(c_double),POINTER(c_double),POINTER(c_double),c_int]

        lib.GalPot_vcsquare.restype = c_double
        lib.GalPot_vcsquare.argtypes = [c_void_p,POINTER(c_double),POINTER(c_double),c_int]

        lib.GalPot_Mass.restype = c_double
        lib.GalPot_Mass.argtypes = [c_void_p,POINTER(c_double),POINTER(c_double),c_int]


        lib.GalPot_RfromLc_single.restype = c_double
        lib.GalPot_RfromLc_single.argtypes = [c_void_p,c_double]

        lib.GalPot_LfromRc_single.restype = c_double
        lib.GalPot_LfromRc_single.argtypes = [c_void_p,c_double]

        lib.GalPot_KapNuOm_single.restype = c_void_p
        lib.GalPot_KapNuOm_single.argtypes = [c_void_p,c_double,POINTER(c_double),POINTER(c_double),POINTER(c_double)]

        # Convert string to usable format and construct the Galaxy potential
        mutable_string = create_string_buffer(str.encode(filename))
        self.obj = lib.GalPot_new(mutable_string)

    def __del__(self):
        '''Removes pointer declared new in the relevant C++ code'''
        return lib.GalPot_delete(self.obj)

    def Potential(self,R,z):
        """Returns potential at R, z

        Parameters:
            R: Float or array of floats.
                Galactocentric radius (cylindrical) in kpc. 
            z: Float or array of floats.
                Galactocentric z in kpc. 
        Returns:
            pot: Float or numpy array.
                Potential in M_solar kpc**2/Myr**2. 
        """
        if len(np.array(R).reshape(-1)) == 1:
            return lib.GalPot_Potential_single(self.obj,c_double(R),c_double(z)) #* self.kpc_Myr_to_km_s**2
        length = len(R)
        assert len(z) == length
        ArrayWithLength = c_double * length
        pot = ArrayWithLength()
        Raux = ArrayWithLength()
        for i,j in enumerate(R):
            Raux[i] = j
        zaux = ArrayWithLength()
        for i,j in enumerate(z):
            zaux[i] = j
        #lib.GalPot_Potential(self.obj,byref(pot),byref(Raux),byref(zaux),c_int(length))
        lib.GalPot_Potential(self.obj,pot,Raux,zaux,c_int(length))
        pot = np.array([p for p in pot])
        return pot #* self.kpc_Myr_to_km_s**2

    def Potential_derivatives(self,R,z):
        """Returns potential and its derivatives wrt R,z at R, z

        Parameters:
            R: Float or array of floats.
                Galactocentric radius (cylindrical) in kpc. 
            z: Float or array of floats.
                Galactocentric z in kpc. 

        Returns:
            pot: Float or numpy array.
                Potential in M_solar kpc**2/Myr**2. 
            dPdR: Float or numpy array.
                Derivative of pot wrt to R in M_solar kpc/Myr**2. 
            dPdz: Float or numpy array.
                Derivative of pot wrt to z in M_solar kpc/Myr**2. 
        """
        if len(np.array(R).reshape(-1)) == 1:
            dPdR = c_double()
            dPdz = c_double()
            pot = lib.GalPot_Potential_and_derivatives_single(self.obj,c_double(R),c_double(z),dPdR,dPdz)
            pot = np.float64(pot)
            dPdR = np.float64(dPdR)
            dPdz = np.float64(dPdz)
            return pot , dPdR, dPdz
            #return pot * self.kpc_Myr_to_km_s**2, dPdR * self.kpc_Myr_to_km_s**2, dPdz * self.kpc_Myr_to_km_s**2

        length = len(R)
        assert len(z) == length
        ArrayWithLength = c_double*length
        pot = ArrayWithLength()
        dPdR = ArrayWithLength()
        dPdz = ArrayWithLength()
        Raux = ArrayWithLength()
        for i,j in enumerate(R):
            Raux[i] = j
        zaux = ArrayWithLength()
        for i,j in enumerate(z):
            zaux[i] = j
        #lib.GalPot_Potential_and_derivatives(self.obj,byref(pot),byref(Raux),byref(zaux),\
        #                                     byref(dPdR),byref(dPdz),c_int(length))
        lib.GalPot_Potential_and_derivatives(self.obj,pot,Raux,zaux,\
                                             dPdR,dPdz,c_int(length))
        pot =np.array([p for p in pot]) #* self.kpc_Myr_to_km_s**2
        dPdR =np.array([p for p in dPdR]) #* self.kpc_Myr_to_km_s**2
        dPdz =np.array([p for p in dPdz]) #* self.kpc_Myr_to_km_s**2
        return pot,dPdR,dPdz

    def ForceRz(self,R,z) :
        """Returns force in R,z directions from potential at R,z

        Force_x = -dP/dx

        Parameters:
            R: Galactocentric radius (cylindrical) in kpc. Float or array of floats.
            z: Galactocentric z in kpc. Float or array of floats.

        Returns:
            Force_R: Force in R direction in M_solar kpc/Myr**2. Float or numpy array.
            Force_z: Force in z direction in M_solar kpc/Myr**2. Float or numpy array.
        """
        P,dPdR,dPdz = self.Potential_derivatives(R,z)
        return -dPdR,-dPdz

    def Acceleration(self,xv_cyl):
        """Returns acceleration in cylindrical coordinates of phase space point xv_cyl in potential

        Note that units are kpc and Myr (and therefore kpc/Myr and kpc/Myr**2)
        Parameters:
            xv_cyl: Coordinates in cylindrical polar coordinates (R,z,phi, v_R,v_z,v_phi)

        Returns:
            acc: acceleration (acc_R, acc_z, acc_phi)
        """
        assert len(xv_cyl) == 6, "xv_cyl not 6 values in acceleration()"
        #order: R, z, phi, vR, vz, vphi
        force = self.ForceRz(xv_cyl[0],xv_cyl[1])
        acc = np.array([force[0]+xv_cyl[5]**2/xv_cyl[0],
                        force[1],
                        -xv_cyl[3]*xv_cyl[5]/(xv_cyl[0])])
        acc = np.array([force[0]+xv_cyl[5]**2/xv_cyl[0],
                        force[1],
                        -xv_cyl[3]*xv_cyl[5]/(xv_cyl[0])])
        return acc

    def XVDerivative(self,time,xv_cyl):
        """Returns time derivatives in cylindrical coordinates of phase space point xv_cyl in potential

        Note that units are kpc and Myr (and therefore kpc/Myr and kpc/Myr**2)
        This function is convenient when integrating an orbit
        Parameters:
            t: Time. Ignored.
            xv_cyl: Coordinates in cylindrical polar coordinates (R,z,phi, v_R,v_z,v_phi)

        Returns:
            deriv: d(xv_cyl)/dt, in cylindrical polar coordinates d(R,z,phi, v_R,v_z,v_phi)/dt
        """
        assert len(xv_cyl) == 6, "xv_cyl not 6 values in acceleration()"
        #order: R, z, phi, vR, vz, vphi
        deriv = np.zeros(6)
        deriv[0:2] = xv_cyl[3:5]
        deriv[2]  = xv_cyl[5]/xv_cyl[0]
        deriv[3:] = self.Acceleration(xv_cyl)[:]
        return deriv

    def Density(self,R,z):
        """Returns total density at R, z

        Parameters:
            R: Galactocentric radius (cylindrical) in kpc. Float or array of floats.
            z: Galactocentric z in kpc. Float or array of floats.

        Returns:
            rho: density in  M_solar/kpc**3 Float or numpy array.

        """
        if len(np.array(R).reshape(-1)) == 1:
            return lib.GalPot_Density_single(self.obj,c_double(R),c_double(z))
        length = len(R)
        assert len(z) == length
        ArrayWithLength = c_double*length
        rho = ArrayWithLength()
        Raux = ArrayWithLength()
        for i,j in enumerate(R):
            Raux[i] = j
        zaux = ArrayWithLength()
        for i,j in enumerate(z):
            zaux[i] = j
        lib.GalPot_Density(self.obj,rho,Raux,zaux,c_int(length))
        #lib.GalPot_Density(self.obj,byref(rho),byref(Raux),byref(zaux),c_int(length))
        rho = np.array([p for p in rho])
        return rho

    def Vcirc2(self,R):
        """Returns velocity squared of circular orbit of radius R

        Parameters:
            R: Galactocentric radius (cylindrical) in kpc. Float or array of floats.

        Returns:
            vc: velocity squared of circular orbit in kpc**2/Myr**2. Float or numpy array.

        """
        if len(np.array(R).reshape(-1)) == 1:
            return lib.GalPot_vcsquare_single(self.obj,c_double(R)) #*self.kpc_Myr_to_km_s**2
        length = len(R)
        ArrayWithLength = c_double*length
        vc2 = ArrayWithLength()
        Raux = ArrayWithLength()
        for i,j in enumerate(R):
            Raux[i] = j
        lib.GalPot_vcsquare(self.obj,vc2,Raux,c_int(length))
        #lib.GalPot_vcsquare(self.obj,byref(vc2),byref(Raux),c_int(length))
        vc2 = np.array([p for p in vc2])
        return vc2 #* self.kpc_Myr_to_km_s**2

    def Vcirc(self,R):
        """Returns velocity of circular orbit of radius R

        Parameters:
            R: Galactocentric radius (cylindrical) in kpc. Float or array of floats.

        Returns:
            vc: velocity of circular orbit in kpc/Myr. Float or numpy array.
        """
        return np.sqrt(self.Vcirc2(R))

    def Mass(self,R):
        """Returns Mass interior to radius R (if 0, total mass)

        Parameters:
            R: Galactocentric radius (cylindrical) in kpc. Float or array of floats.

        Returns:
            mass: mass interior to R in M_solar. Float or numpy array.

        """
        if len(np.array(R).reshape(-1)) == 1:
            return lib.GalPot_Mass_single(self.obj,c_double(R))
        length = len(R)
        ArrayWithLength = c_double*length
        mass = ArrayWithLength()
        Raux = ArrayWithLength()
        for i,j in enumerate(R):
            Raux[i] = j
        lib.GalPot_Mass(self.obj,mass,Raux,c_int(length))
        #lib.GalPot_Mass(self.obj,byref(mass),byref(Raux),c_int(length))
        mass = np.array([p for p in mass])
        return mass

    def KapNuOm(self,R):
        """Returns epicyclic frequencies of orbit at radius R

        Parameters:
            R: Galactocentric radius (cylindrical) in kpc. Float.

        Returns:
            kappa: Radial epicycle frequency at R in 1/Myr. Float.
            nu: Vertical epicycle frequency at R in 1/Myr. Float.
            Om: Circular frequency at R in 1/Myr. Float.
        """
        kappa,nu,Om = c_double(0.),c_double(0.),c_double(0.)
        lib.GalPot_KapNuOm_single(self.obj,c_double(R),kappa,nu,Om)
        
        return np.double(kappa),np.double(nu),np.double(Om)

    def RfromLc(self,Lz):
        """Returns Radius of circular orbit with angular momentum Lz

        Parameters:
            Lz: Angular momentum in kpc**2/Myr. Float.

        Returns:
            R: Galactocentric radius (cylindrical) in kpc. Float.
        """
        return np.double(lib.GalPot_RfromLc_single(self.obj,c_double(Lz)))

    def LfromRc(self,R):
        """Returns angular momentum of circular orbit of radius R

        Parameters:
            R: Galactocentric radius (cylindrical) in kpc. Float.

        Returns:
            Lz: Angular momentum in kpc**2/Myr. Float.
        """
        return np.double(lib.GalPot_LfromRc_single(self.obj,c_double(R)))


#---------------------------------------------------------------------------
#
#   OrbitIntegrator class
#
#---------------------------------------------------------------------------


class OrbitIntegrator:

    """Interface for orbit integration in GalPot

    Parameters:
    ----------
    GalaxyPotentialInput : GalaxyPotential
        GalaxyPotential in which the orbit will be integrated.
    

    This is put here as a convenient method using the orbit integrator with 
    GalPot in Python.

    The original C++ code is still the main element, and this just calls 
    the relevant functions. It is noticibly faster than using Python 
    integrators if you are only interested in the orbital parameters.

    Everything is given in GalaxyPotential's internal units, with are 
    kpc, Myr and M_solar.
    To convert to km/s, it may be convenient to use the value 
    GalaxyPotential.kpc_Myr_to_km_s 

    Methods
    ----------
    getOrbitStats(XV, t_end) :      returns OrbitStat, a named tuple 
                                    containing orbit properties
    getOrbitPathandStats(XV,times): returns position and velocity on orbit 
                                    at specified times, and OrbitStat

    """

    kpc_Myr_to_km_s = 977.77

    def __init__(self,GalaxyPotentialInput):
        """Intialise interface for orbit integration in GalPot

        Parameters:
            GalaxyPotentialInput : GalaxyPotential
                GalaxyPotential in which the orbit will be integrated.
        """
        
        # Declare various things from the C++ code
        lib.OrbitIntegrator_new.restype = c_void_p
        lib.OrbitIntegrator_new.argtypes = [c_void_p, c_double]
        lib.OrbitIntegrator_delete.restype = c_void_p
        lib.OrbitIntegrator_delete.argtypes = [c_void_p]
        lib.resetOrbitIntegratorEndTime.restype = c_void_p
        lib.resetOrbitIntegratorEndTime.argtypes = [c_void_p, c_double]
        lib.runOrbitIntegratorforstats.restype = c_int
        lib.runOrbitIntegratorforstats.argtypes =[c_void_p,POINTER(c_double),
                                                  POINTER(c_double),c_int]
        lib.runOrbitIntegratorwithPath.respath = c_int
        lib.runOrbitIntegratorwithPath.argtypes = [c_void_p,POINTER(c_double),
                                                   POINTER(c_double),c_int,
                                                   POINTER(c_double),
                                                   POINTER(c_double),c_int]

        # Convert string to usable format and construct the Galaxy potential
        self.obj = lib.OrbitIntegrator_new(GalaxyPotentialInput.obj,c_double(13800.))

    def __del__(self):
        '''Removes pointer declared new in the relevant C++ code'''
        return lib.OrbitIntegrator_delete(self.obj)


    def getOrbitStats(self,XV, t_end=13800.):
        """Returns stats for the orbit that starts at point XV

        Parameters:
            XV:     ndarray or list (dimensions N x 6)
                Position & velocity [R,z,phi,v_R,v_z,v_phi] with units 
                kpc, radians, kpc/Myr.
            t_end: float (default 13800.)
                Final integration time (by default 13800 Myr)

        Returns:

            Energy: ndarray (or float if one dimensional array input)           
                Orbital Energy. If positive, orbit is unbound, and 
                values other than Lz are set to 0
            Lz: ndarray (or float if one dimensional array input)
                Angular Momentum
            Pericentre: ndarray (or float if one dimensional array input)
                Minimum Galactocentric radius (r = sqrt(R**2+z**2))
            Apocentre: ndarray (or float if one dimensional array input)
                Maximum Galactocentric radius (r = sqrt(R**2+z**2))
            Zmax: ndarray (or float if one dimensional array input)
                Maximum Galactocentric height
            GuidingRadius: ndarray (or float if one dimensional array input)
                Radius of a circular orbit with the same angular momentum
            PseudoEccentricity: ndarray (or float if one dimensional array input)
                (Apocentre-Pericentre)/(Apocentre+Pericentre). I call it 
                "Pseudo" becuase eccentricity only really makes sense in a
                sperically symmetric system

        """
        
        oneD = False
        
        try:
        # `sample` is an ND-array.
            Nvalues, Ndim = XV.shape
        except (AttributeError, ValueError):
        # `sample` is (a sequence of) 1D arrays.
            XV = np.atleast_2d(XV)
            Nvalues, Ndim = XV.shape
            if Nvalues == 1 : oneD=True
        
        assert (Ndim == 6), ("Input must have 6 values per array row "
                             +"(R,z,phi,vR,vz,vphi)")

        ArrayWithLength6N = c_double * (6*Nvalues)
        ArrayWithLength7N = c_double * (7*Nvalues)
        XV_c = ArrayWithLength6N()
        output_c = ArrayWithLength7N()
        for i,XVval in enumerate(XV.flatten()):
            XV_c[i] = XVval
        
        #lib.resetOrbitIntegratorEndTime(self.obj,c_double(t_end))
        lib.resetOrbitIntegratorEndTime(self.obj,c_double(t_end))
        lib.runOrbitIntegratorforstats(self.obj,XV_c,output_c,Nvalues)
        
        Energy              = output_c[0::7]
        Lz                  = output_c[1::7]
        Pericentre          = output_c[2::7]
        Apocentre           = output_c[3::7]
        Zmax                = output_c[4::7]
        GuidingRadius       = output_c[5::7]
        PseudoEccentricity  = output_c[6::7]
        OrbitStat = namedtuple('OrbitStat',
                               ('Energy','Lz','Pericentre','Apocentre','Zmax',
                                'GuidingRadius','PseudoEccentricity'))
        if oneD : 
            return OrbitStat(Energy[0],Lz[0],Pericentre[0],Apocentre[0],
                             Zmax[0],GuidingRadius[0],PseudoEccentricity[0])
        return OrbitStat(np.array(Energy),np.array(Lz),np.array(Pericentre),
                         np.array(Apocentre),np.array(Zmax),np.array(GuidingRadius),
                         np.array(PseudoEccentricity))


    def getOrbitPathandStats(self,XV, times):
        """Returns orbital path and stats for the orbit that starts at point XV

        Orbit is integrated over the times given, which should be in 
        strictly increasing (or decreasing) order. Times are in Myr.
        Parameters:
            XV:     ndarray or list (dimensions N x 6)
                Position & velocity [R,z,phi,v_R,v_z,v_phi] with units 
                kpc, radians, kpc/Myr.

            times - ndarray or list (1 dimensional)
                Times for output of orbital path, with units Myr
       
        Returns:
            paths:
                Position & velocity [R,z,phi,v_R,v_z,v_phi] on the orbits at 
                the specified times. 

            OrbitStats: a named tuple containing the following
                Energy: ndarray (or float if one dimensional array input)           
                    Orbital Energy. If positive, orbit is unbound, and 
                    values other than Lz are set to 0
                Lz: ndarray (or float if one dimensional array input)
                    Angular Momentum
                Pericentre: ndarray (or float if one dimensional array input)
                    Minimum Galactocentric radius (r = sqrt(R**2+z**2))
                Apocentre: ndarray (or float if one dimensional array input)
                    Maximum Galactocentric radius (r = sqrt(R**2+z**2))
                Zmax: ndarray (or float if one dimensional array input)
                    Maximum Galactocentric height
                GuidingRadius: ndarray (or float if one dimensional array input)
                    Radius of a circular orbit with the same angular momentum
                PseudoEccentricity: ndarray (or float if one dimensional array input)
                    (Apocentre-Pericentre)/(Apocentre+Pericentre). I call it 
                    "Pseudo" becuase eccentricity only really makes sense in a
                    sperically symmetric system

        """
        # find maximum (or minimum if negative)
        t_end = times[np.argmax(np.absolute(times))]

        oneD = False
        
        try:
        # `sample` is an ND-array.
            Nvalues, Ndim = XV.shape
        except (AttributeError, ValueError):
        # `sample` is (a sequence of) 1D arrays.
            XV = np.atleast_2d(XV)
            Nvalues, Ndim = XV.shape
            if Nvalues == 1 : oneD=True
        nt = len(times)
        assert (Ndim == 6), ("Input must have 6 values per array row "
                             +"(R,z,phi,vR,vz,vphi)")

        # Needed for c code
        ArrayWithLength6N = c_double * (6*Nvalues)
        ArrayWithLength7N = c_double * (7*Nvalues)
        ArrayWithLengthnt = c_double * (nt)
        ArrayWithLengthntN6 = c_double * (nt*Nvalues*6)
        XV_c = ArrayWithLength6N()
        input_t_c = ArrayWithLengthnt()
        output_c = ArrayWithLength7N()
        output_path_c = ArrayWithLengthntN6()

        for i,t in enumerate(times) :
            input_t_c[i] = t 
        for i,XVval in enumerate(XV.flatten()):
            XV_c[i] = XVval
        
        # call c++ routines

        lib.resetOrbitIntegratorEndTime(self.obj,c_double(t_end))

        lib.runOrbitIntegratorwithPath(self.obj,XV_c,input_t_c,nt, 
                                       output_path_c,output_c, Nvalues)
        # organise output
        paths = np.array(output_path_c).reshape([Nvalues,-1,6])
        
        Energy              = output_c[0::7]
        Lz                  = output_c[1::7]
        Pericentre          = output_c[2::7]
        Apocentre           = output_c[3::7]
        Zmax                = output_c[4::7]
        GuidingRadius       = output_c[5::7]
        PseudoEccentricity  = output_c[6::7]       
        OrbitStat = namedtuple('OrbitStat',
                               ('Energy','Lz','Pericentre','Apocentre','Zmax',
                                'GuidingRadius','PseudoEccentricity'))
        if oneD : 
            return paths[0],OrbitStat(Energy[0],Lz[0],Pericentre[0],Apocentre[0],
                             Zmax[0],GuidingRadius[0],PseudoEccentricity[0])
        return paths, OrbitStat(np.array(Energy),np.array(Lz),np.array(Pericentre),
                         np.array(Apocentre),np.array(Zmax),np.array(GuidingRadius),
                         np.array(PseudoEccentricity))



