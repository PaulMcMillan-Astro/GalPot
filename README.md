GalPot
============

This is a stand-alone version of [Walter Dehnen's GalaxyPotential C++ code](http://ukads.nottingham.ac.uk/abs/1998MNRAS.294..429D). It is a version I took from the falcON code in the [NEMO Stellar Dynamics Toolbox](http://chara.astro.umd.edu/nemo/), and cut down to stand alone. It is a convenient way of finding the gravitational potential associated with axisymmetric density profiles.

It should compile if you run make.

The potentials provided in pot/ are primarily from McMillan (2016), but also from Piffl et al (2014), McMillan (2011) and Dehnen & Binney (1998).

It is under the GNU General Public Licence (GPL).

It comes with no warranty, especially not from Walter. I leave it here simply for convenience.

The executable file testGalPot.exe shows examples of using GalPot

The other executable files provide simple tools for integration of orbits in the given potentials. Running the cshell script Example.csh (e.g. "csh Example.csh") runs example versions of these executables.

Documentation.pdf provides an explanation of the input files used.

Running any of the executables without any arguments should provide a description of how it is used.

If you use this code, please cite whoever determined the potential (preferably McMillan 2016!), and Dehnen & Binney (1998).