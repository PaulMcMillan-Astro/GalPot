

from ctypes import *
import numpy as np
import os.path

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

    This is put here as a (hopefully) convenient method of using GalPot in this new Pythonic world.
    The original C++ code is still the main element.

    Everything is given in GalaxyPotential's internal units, with are kpc, Mpc and M_solar

    Methods
    ----------
    Potential(R,z) :            returns P, Potential at postion(s)
    Potential_derivatives(R,z): returns P, dPdR, dPdz: Potential and derivatives at position(s)
    ForceRz(R,z) :              returns f_R,f_z, specific force from potential at position(s)
    Acceleration(xv_cyl) :      returns a_R, a_z, a_phi: acceleration of moving body in potential
    XVDerivative(xv_cyl):       returns dxv/dt, derivative of phase space position (required for orbit integration)
    Density(R,z):               returns rho, density at position(s)
    Vcirc2(R):                  returns vcirc^2: square of circular velocity at R (single value or array)
    Vcirc(R):                   returns vcirc: circular velocity at R (single value or array)
    KapNuOm(R)                  returns kappa,nu,omega: epicycle frequencies at R (single value only)
    RfromLc(Lz)                 returns R, radius of circular orbit with angular momentum Lz
    LfromRc(R)                  returns Lz, angular momentum of circular orbit with radius R

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
                                                                POINTER(c_double),POINTER(c_double),c_int]

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
            R: Galactocentric radius (cylindrical) in kpc. Float or array of floats.
            z: Galactocentric z in kpc. Float or array of floats.

        Returns:
            pot: Potential in M_solar kpc**2/Myr**2. Float or numpy array.
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
            R: Galactocentric radius (cylindrical) in kpc. Float or array of floats.
            z: Galactocentric z in kpc. Float or array of floats.

        Returns:
            pot: Potential in M_solar kpc**2/Myr**2. Float or numpy array.
            dPdR: Derivative of pot wrt to R in M_solar kpc/Myr**2. Float or numpy array.
            dPdz: Derivative of pot wrt to z in M_solar kpc/Myr**2. Float or numpy array.
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
        #lib.GalPot_KapNuOm_single(self.obj,c_double(R),byref(kappa),byref(nu),byref(Om))
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
