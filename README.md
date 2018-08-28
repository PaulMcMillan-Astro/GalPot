GalPot
============

This is a stand-alone version of [Walter Dehnen's GalaxyPotential C++ code](http://ukads.nottingham.ac.uk/abs/1998MNRAS.294..429D), with an added Python wrapper for those who want it (see below for details). It is a version I took from the falcON code in the [NEMO Stellar Dynamics Toolbox](http://chara.astro.umd.edu/nemo/), and cut down to stand alone. It is a convenient way of finding the gravitational potential associated with axisymmetric density profiles.

The package also includes code that performs transformations between commonly used coordinate systems for both positions and velocities (the class OmniCoords), and that integrates orbits in the potentials.

It should compile (including the library required for the Python wrapper) if you run make.

The potentials provided in pot/ are primarily from McMillan (2017), but also from Piffl et al (2014), McMillan (2011) and Dehnen & Binney (1998). See pot/README for details. Note that the naming convention for the McMillan (2017) potentials is slightly muddled because the paper and code were originally made available in 2016 (not 2017).

It is under the GNU General Public Licence (GPL).

It comes with no warranty, especially not from Walter. I leave it here simply for convenience.

The executable file testGalPot.exe shows examples of using GalPot

The other executable files provide simple tools for integration of orbits in the given potentials. Running the cshell script Example.csh (e.g. "csh Example.csh") runs example versions of these executables.

In general: Those called findOrbit*.exe require the starting points given on the command line, while those called findManyOrbit*.exe require starting points listed in a file.

Documentation.pdf provides an explanation of the input files used.

Running any of the executables without any arguments should provide a description of how it is used.

If you use this code, please cite whoever determined the potential (preferably [McMillan 2017](http://dx.doi.org/10.1093/mnras/stw2759)!), and [Dehnen & Binney (1998)](http://dx.doi.org/10.1046/j.1365-8711.1998.01282.x).


# New

I have now added the possibility of adding a Kepler potential or Miyamoto-Nagai potential to the standard ones produced by GalPot. The class MultiPotential exists to combine these, and the file findOrbitMultiPot.cc exists to demonstrate this. The other routines have not been updated to allow this functionality, and findOrbitMultiPot should be considered a template which the interested user can follow.

Note that the basic GalPot is rather poorly suited to representing smaller-scale objects (like Sgr A* or a nuclear disc). It is recommended to combine GalPot with a Kepler or Miyamoto-Nagai potential if you want to represent these objects.

## New 8/2018:

I have now added a Python wrapper - note that only the gravitational potential code has been wrapped, because astropy & scipy already provide convenient methods for the other activities.

The library (obj/libPyGalPot.so) is compiled with make, along with the C++ executables. The wrapper is then in the file GalPot.py

GalPotPythonExample.ipynb demonstrates the Python wrapper. The orbit integration is done with scipy.integrate.solve_ivp, which is new in scipy 1.0.0 - that part will only work if your scipy install is up to date.

The intention of the Python wrapper is that it is available for people who want to do things that are more complicated than the things already done in the executables, but don't want to work in c++.

Warning: it is significantly slower than the raw c++ code.

Comments on usability welcome.
