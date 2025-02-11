# Read me for collion.cpp

A gas in collisional ionization equilibrium

This program varies the gas kinetic temperature in an environment where collisional ionization will dominate.  This is very similar to conditions in the solar corona. 

The gas temperature ranges between 0.5 < log T < 9.5 in 0.1 dex steps.  This covers the code's intended temperature range. The hydrogen density is 1.0 cm-3.  
Cosmic rays are added so that the ion-molecule UMIST chemistry can function. We want to test thermal processes so the CR density set to a small value. 
This is needed for the chemistry to function.

The predicted ionization fractions are placed in the file collion.txt which is formatted to be similar to Jordan 1969, MNRAS, 142, 501.
The script split-collion.pl will spit collion.txt into separate files for each element.

one.in has the program commands in an input deck that can be run stand-alone for testing.
