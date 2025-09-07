## Notes

This directory contains the files needed to generate
a clang++ executable that is optimized for line-level profiling

Edit Makefile.conf so that the default build has other properties,
or use this to pass additional options to make:
```
make -j <n> EXTRA="-DFOO -DBAR"
```
