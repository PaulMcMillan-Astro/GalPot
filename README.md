GalPot
============

This is a stand-alone version of [Walter Dehnen's GalaxyPotential C++ code](https://ui.adsabs.harvard.edu/abs/1998MNRAS.294..429D), with an added Python wrapper for those who want it (see below for details). It is a version I took from the falcON code in the [NEMO Stellar Dynamics Toolbox](https://github.com/teuben/nemo), and cut down to stand alone. It is a convenient way of finding the gravitational potential associated with axisymmetric density profiles. It now comes with a Python wrapper.

The package also includes code that performs transformations between commonly used coordinate systems for both positions and velocities (the class OmniCoords), that integrates orbits in the potentials, and that provides other potentials (Kepler potential, Miyamoto-Nagai potential, or a combination of any number of these and GalaxyPotential using the class MultiPot).

There is a Python wrapper for the gravitational potential and orbit integration code. Astropy already provides convenient methods for coordinate transforms. The orbit integration is dramatically faster (around 20 times faster in tests) than using a Python integrator.

It should all compile (including the library required for the Python wrapper) if you run make.

The potentials provided in pot/ are primarily from McMillan (2017), but also from Piffl et al (2014), McMillan (2011) and Dehnen & Binney (1998). See pot/README for details. Note that the naming convention for the McMillan (2017) potentials is slightly muddled because the paper and code were originally made available in 2016 (not 2017).

It is under the GNU General Public Licence (GPL).

It comes with no warranty, especially not from Walter. I leave it here simply for convenience.

If you use this code, please cite whoever determined the potential (preferably [McMillan 2017](http://dx.doi.org/10.1093/mnras/stw2759)!), and [Dehnen & Binney (1998)](http://dx.doi.org/10.1046/j.1365-8711.1998.01282.x).

## MCMC chains from McMillan (2017)

The MCMC chain used to explore the probability distribution of the model parameters for the main results given in McMillan (2017) are available [via Zenodo](https://doi.org/10.5281/zenodo.6219842). The notebook ReadAndSelectMCMCOutput.ipynb gives an example of reading the chain and sampling a model from it.

# Examples and documentation

The executable file testGalPot.exe shows examples of using GalPot

The other executable files provide simple tools for integration of orbits in the given potentials. Running the cshell script Example.csh (e.g. "csh Example.csh") runs example versions of these executables.

In general: Those called findOrbit*.exe require the starting points given on the command line, while those called findManyOrbit*.exe require starting points listed in a file.

Documentation.pdf provides an explanation of the input files used.

Running any of the executables without any arguments should provide a description of how it is used.

GalPotPythonExample.ipynb demonstrates the Python wrapper. 


### Notes: 

The basic GalPot is rather poorly suited to representing smaller-scale objects (like Sgr A* or a nuclear disc). It is recommended to combine GalPot with a Kepler or Miyamoto-Nagai potential if you want to represent these objects.


