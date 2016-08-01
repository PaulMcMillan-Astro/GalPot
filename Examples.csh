#!/bin/csh

echo Example script to show the executables available. Run after you have run make

echo See also testGalPot

# Demonstrate findOrbitProperties.exe
echo
echo Find orbital properties of the Sun...
echo Running "./findOrbitProperties.exe 8.2 0.014 11.1 7.25 -245"

./findOrbitProperties.exe 8.2 0.014 11.1 7.25 -245

echo
echo "Finding orbital path of the Sun (5000 points), and putting it in a file..."
echo Running "./findOrbit.exe 8.2 0.014 0 11.1 7.25 -245 examples/Output_SunOrbit.tab 5000"
./findOrbit.exe 8.2 0.014 0 11.1 7.25 -245 examples/Output_SunOrbit.tab 5000

echo
echo "Finding the properties of many orbits from initial conditions specified in ExampleManyOrbitProperties.tab"
echo See examples/Output_ManyOrbitExample.tab for the output...
echo Running "./findManyOrbitProperties.exe examples/ExampleManyOrbitProperties.tab examples/Output_ManyOrbitExample.tab"
./findManyOrbitProperties.exe examples/ExampleManyOrbitProperties.tab examples/Output_ManyOrbitExample.tab

echo
echo "Finding properties of orbits assuming they are in Equatorial coordinates"
echo See examples/Output_FromEqua.tab for output
echo Running "./findManyOrbitPropertiesfromEquatorial.exe examples/ExampleManyFromObservable.tab examples/Output_FromEquat.tab"
./findManyOrbitPropertiesfromEquatorial.exe examples/ExampleManyFromObservable.tab examples/Output_FromEquat.tab

echo
echo "Finding properties of orbits assuming they are in Galactic coordinates"
echo See examples/Output_FromGal.tab for output
echo Running "./findManyOrbitPropertiesfromGalactic.exe examples/ExampleManyFromObservable.tab examples/Output_FromGal.tab"
./findManyOrbitPropertiesfromGalactic.exe examples/ExampleManyFromObservable.tab examples/Output_FromGaltab

echo
echo "Find properties of orbits in Equatorial coordinates, with uncertainties, using 1000 Monte Carlo samples"
echo See examples/Output_WithErrors.tab for output
echo Running "./findManyOrbitPropertiesfromEquatorialwErrors.exe examples/ExampleEquatorialwErrors.tab 1000 examples/Output_WithErrors.tab"
./findManyOrbitPropertiesfromEquatorialwErrors.exe examples/ExampleEquatorialwErrors.tab 1000 examples/Output_WithErrors.tab

echo

#end
