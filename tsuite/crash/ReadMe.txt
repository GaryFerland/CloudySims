These test various abort modes.   The tests are:

assert.in the assert macro

bounds_auto.in test array bounds
bounds_heap.in
bounds_multi.in
bounds_multi_iter.in
bounds_static.in

exception.in throws a C++ exception

fpe_isnan.in floating point error (FPE) is NaN
fpe_longoverflow.in floating point error
fpe_NaN.in floating point error
fpe_overflow.in floating point error
fpe_setnan.in floating point error
fpe_zero.in floating point error

grid_corners.in corners of grid exceed temperature limits of code,
test that control passes smoothly back to main with all punch produced after
abort is declared

undef.in variable is undefined

runall.pl - perlo script to run all tests
