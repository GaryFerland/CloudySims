two tests that crash with different errors
cmesh fails with the energy mesh for x-ray 2-electron

ca1 turns of elements causing above, and so gets to the failure
with Ca I

==========================================

cmesh crashes with
           Tbr(FarIR):1.596E-07   Tbr(H n=6):9.778E-06   Tbr(1Ryd):4.616E-10   Tbr(4R):9.505E-12     Tbr (Nu-hi):1.622E-33
  
SanityCheck found non-positive photo cs for nelem=28 ion=27 
value was 0.00e+00 + 0.00e+00 nelem 28 ion 27 at energy 8.13e+02

-----------------------

ca1.in

 PROBLEM atom_levelN found negative population, nNegPop=7, atom=Ca 1 lgSearch=F Te=  8.71e+03 fnzone 0.00 
 Absolute:   1.29e-03  9.53e-08  2.78e-07  4.38e-07 -1.74e-06 -2.89e-06 -4.00e-06  2.01e-05  3.09e-09  1.10e-11  4.49e-07 -1.89e-10 -2.53e-10 -1.79e-10 -3.12e-10
 Relative:   9.90e-01  7.31e-05  2.13e-04  3.36e-04 -1.34e-03 -2.21e-03 -3.07e-03  1.54e-02  2.37e-06  8.44e-09  3.44e-04 -1.45e-07 -1.94e-07 -1.37e-07 -2.39e-07
 PROBLEM in dBase_solve, atom_levelN returned negative population .
 No parts of the continuum were negative, the electron density was  2.74e+03 te=  8.71e+03
 This is zone number   0

