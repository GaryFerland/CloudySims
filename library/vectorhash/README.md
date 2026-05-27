## VectorHash - a very fast vectorized checksum algorithm

VectorHash is an extremely fast, non-cryptographic checksum algorithm. It is
primarily designed to create high quality checksums of files or buffers, but can
be used for other purposes as well. In this project the term hash will be used
as a synonym for checksum, even though the term hash has a broader definition
than checksum. In principle this algorithm can also be used to populate hash
tables, but in practice it will likely be too slow for that purpose if the keys
are small in size. This is because the minimum blocksize that the code can work
with is rather large compared to other algorithms.

This program has been written such that it can be vectorized using SIMD
instructions from the SSE2, AVX2, or AVX512f instruction sets, which makes it
very fast. It also has a scalar version for non-Intel platforms. The resulting
checksum is identical for each of these instruction sets. The fact that SIMD
instructions work on very wide registers, combined with the requirement that all
versions of the algorithm create identical checksums, results in the fact that
the minimum blocksize that the algorithm can handle is rather large. So the very
thing that makes the algorithm extremely fast on larger buffers, creates a lot
of overhead for very small keys. The algorithm will show good speed for buffers
roughly 16--32 kiB and larger. A more detailed breakdown of the speed of the
algorithm can be seen below.

<div style="text-align: center;">
<img src="./speedtest.png" width="50%"/>
</div>

This image shows the bulk speed of the algorithm (in GiB/s) as a function of the
buffer size. The latter is shown as the base-2 logarithm of the size in kiB, so
the left-most point in the plot is for a buffer of 2<sup>0</sup> = 1 kiB and the
right-most point for a buffer of 2<sup>12</sup> = 4096 kiB = 4 MiB. Curves are
shown for all 4 versions of the algorithm: scalar, SSE2, AVX2, and AVX512f. It
can be seen that the scalar and SSE2 versions have nearly identical speeds,
while the AVX2 and AVX512f versions give considerable speedups. For the AVX2 and
AVX512f versions, the speed peaks for a buffer size of 512 kiB and then drops
substantially for larger buffers. Most likely this is caused by the size of the
CPU cache: for very large buffers the speed is limited by direct memory access,
while for smaller buffers it is limited by the speed of the CPU cache. These
tests were carried out with release v1.1 on an Intel(R) Xeon(R) Gold 6226R CPU
@ 2.90GHz.

VectorHash has in part been inspired by the PRNG algorithm
[xoroshiro128++](https://prng.di.unimi.it/) written by David Blackman and
Sebastiano Vigna (for the core of the hashing function) and
[MurmurHash3](https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp)
written by Austin Appleby (for the finalization mix and some support routines).

The algorithm can produce checksums for many different widths (any multiple of
32 bits between 32 and 1024 bits). It has been extensively tested against
SMHasher written by ZhuReini Urban for all the widths listed below. The code
passes all the tests using commit 83fd13bc of SMHasher, except for the speed
test. The latter fails because it measures the speed for hashing very small
keys, which is not what the algorithm is designed for as was already explained
above. The test harness can be found [in this repository](https://github.com/rurban/smhasher).
The output of the test runs is here:
[32bit](https://gitlab-as.oma.be/-/project/876/uploads/df1fdf1da4d3d89fa6b1e30ea581d831/test_results_v1.0_32bit.txt),
[64bit](https://gitlab-as.oma.be/-/project/876/uploads/51c25fc5b5b05480a5324bcb8de857c6/test_results_v1.0_64bit.txt),
[96bit](https://gitlab-as.oma.be/-/project/876/uploads/2b09dffff73ff07217703b2a41dc0a5c/test_results_v1.0_96bit.txt),
[128bit](https://gitlab-as.oma.be/-/project/876/uploads/656fc6555fffcfa56f0592fbadd64f8d/test_results_v1.0_128bit.txt),
[160bit](https://gitlab-as.oma.be/-/project/876/uploads/bc4a089acde079baa854a234f0d12354/test_results_v1.0_160bit.txt),
[192bit](https://gitlab-as.oma.be/-/project/876/uploads/2caf1884915025aa5abc6bd07ec121a0/test_results_v1.0_192bit.txt),
[224bit](https://gitlab-as.oma.be/-/project/876/uploads/3b43648216037ff5f1896e04445b09c6/test_results_v1.0_224bit.txt),
[256bit](https://gitlab-as.oma.be/-/project/876/uploads/ecf2c66de7faeeb077291c937d19bb4f/test_results_v1.0_256bit.txt),
[512bit](https://gitlab-as.oma.be/-/project/876/uploads/af5609c5ac414e3bcdc7a23a350ac5da/test_results_v1.0_512bit.txt),
[1024bit](https://gitlab-as.oma.be/-/project/876/uploads/60d84a7037cfb5b6f12c41ef4577db8d/test_results_v1.0_1024bit.txt).

### Command line routines

After building and installing the code (see the separate file INSTALL.md for
this) you will have three executables called vh128sum, vh256sum, and vh512sum.
These can produce checksums of any given file. They are designed to be a plug-in
replacement for other well-known checksum algorithms. The actual checksums will
be different of course, but the aim is to replicate all command line options.
This makes it easy to adapt scripts, etc. See the man page included in this
distribution for a more detailed description of the options. In the remaining
description only the name vh128sum will be used for brevity, but all remarks
apply implicitly to the other executables as well.

The executables are in fact identical and are hard links to each other. The code
determines the width of the checksum by looking at the name of the executable.
If that name is of the form "vh\<number\>sum" it will produce a \<number\>-bit
wide checksum, e.g. an executable named vh256sum will produce a 256-bit wide
checksum. This works for all supported checksum widths. If it is not possible to
determine the width of the checksum this way, the width will default to 128
bits. It is also possible to force the width with command line flags, e.g. the
command vh128sum \-l256 will produce a 256-bit checksum. Supported lengths are
any multiple of 32 bits between 32 and 1024. Not all checksum widths have been
tested against SMHasher. See the previous section for a list of all the tests.

The vh128sum executable can handle multiple file names in a single invocation
(as long as you do not exceed the maximum command-line length of the shell).
Creating or checking checksums this way will be much faster than invoking
vh128sum separately for each individual file. So this practice is strongly
recommended.

The code will automatically detect if supported SIMD instructions are available
and use the highest set available. If you include the \--verbose flag, the code
will report which SIMD instruction set has been used. You can also override this
detection with command line options, e.g. vh128sum \--avx512 will force the code
to use the AVX512f instruction set, even if the hardware does not support those
instructions. You will get a crash on an illegal instruction in the latter case.
These flags are intended for testing and debugging and should not be needed in a
normal production environment.

The command vh128sum \--help will give a complete overview of all flags that are
supported.

The checksums of an empty file are as follows. These can be reproduced with the
command

	echo -n | vh128sum -l 32

and similarly for other checksum widths.

32 bit:  
<tt>7647d9bd</tt>

64 bit:  
<tt>73711a77d6031b6f</tt>

96 bit:   
<tt>fdcdf538d6031b6f8a613d7f</tt>

128 bit:   
<tt>fe82e7d9998e9819c7ac954ea0a0ea8e</tt>

160 bit:   
<tt>a49af44fdf9952d380d2cb84aace19bf6c0f58c8</tt>

192 bit:   
<tt>5481034fd2a9074478684acf85dceecd94e0e7830f82775b</tt>

224 bit:   
<tt>cfd4bf9386cf00766add0e4a30fc6e7ca6cf10f922cdf8e6cc5578ac</tt>

256 bit:   
<tt>c87638f19f1778d7ee03dddc5c6677fd0701e532d311bfffddfc21c0ece8910a</tt>

512 bit:   
<tt>280f647d93f843b75627cc8c437084e2668a3cd9e077a2cf3186c67c990983c2</tt>
<tt>8ff8180187d362a79cb02cb49ba350a5e185e0050fc273e836b1db30571cfa66</tt>

1024 bit:   
<tt>c38bd292e55400931e47f4397b5266b3e143ffd918fb944c33ad908458a08c62</tt>
<tt>65e9c28007bf68301ab6f68c82ddcc2fa993ad69556b9ab4d588594ac4d02b5d</tt>
<tt>dd083a4cf1418b6522f06bff90a67bf93faaa719e1dbdbda3592deddbd1a8378</tt>
<tt>f0993653e73c292f338729c1efe01b352e90ed9a0ae01e613ad72736db619b2d</tt>

### The libvhsum library

After building and installing the code you will have a 64-bit (and/or 32-bit)
shared library called <tt>libvhsum.so</tt>. This is a symbolic link to a
versioned library called <tt>libvhsum.so.n</tt> where <tt>n</tt> is the major
version number of the release (e.g., <tt>libvhsum.so.1</tt>). Different releases
of the library with the same major version number will be backward compatible.
This library exports two routines that will allow you to perform a checksum on a
buffer. The declaration for the call is:

    #include <vectorhash.h>
    void VectorHash(const void *buf, size_t len, uint32_t seed, void *out, size_t hw);
    void VectorHashSIMD(const void *buf, size_t len, uint32_t seed, void *out, is_type SIMDtype, size_t hw);

Here the variable <tt>buf</tt> is a pointer to the buffer that should be
checksummed, <tt>len</tt> is the length of the buffer, <tt>seed</tt> is the seed
for the checksum algorithm (use 0xfd4c799d to replicate the behavior of
vh128sum, etc, but any other value is fine too), <tt>out</tt> is a pointer to
the buffer that will be used to write the checksum into, and <tt>hw</tt> is the
width of the checksum in bits (allowed values are any multiple of 32 between 32
and 1024). The routine VectorHashSIMD() allows you to specify the SIMD
instruction set to be used (allowed values are: IS_SCALAR, IS_SSE2, IS_AVX2, and
IS_AVX512). This version should not normally be used, but can be helpful for
e.g. speed tests. The declaration is contained in the header file
<tt>vectorhash.h</tt>. A very simple program using this library could contain:

    #include <vectorhash.h>
    // ...
    size_t len = 1048576; // 1 MiB, but any other size is fine too
    void *buf = malloc(len);
    // fill buffer...
    static const size_t hw = 256; // get a 256-bit hash
    uint32_t checksum[hw/32];
    VectorHash(buf, len, 0xfd4c799d, checksum, hw);

and would need to be compiled with a command like this:

    g++ test.cc -lvhsum -L/path/to/libdir

The <tt>-L</tt> parameter is only needed if the library is not installed in a
directory that is included in your <tt>LIBRARY_PATH</tt>. The library is
written in C++, but can also be called from C programs (and even other languages
like e.g. Fortran) using one additional parameter:

	gcc test.c -lvhsum -L/path/to/libdir -lstdc++

The routine <tt>VectorHash</tt> will automatically determine the hardware
capabilities of the processor and use the appropriate version of the algorithm.
However, there is one caveat. Assuming your processor supports AVX512f
instructions, it will need to read 512 bits of data at once from memory. For
this the memory address (and hence also the start of the buffer) will need to be
aligned on a 512-bit = 64-byte boundary (in other words, the last 6 bits of the
buffer address need to be zero). For AVX2 instructions 32-byte alignment is
needed, and for SSE2 instructions 16-byte alignment. The scalar version can work
on any buffer alignment (even odd addresses). If the buffer is not correctly
aligned, the code will automatically revert to a lower algorithm that does
support the provided buffer alignment. For instance, if the hardware supports
AVX512f instructions, but the buffer is aligned on a 16-byte boundary, the code
will actually use the SSE2 version of the algorithm since using the AVX512f or
AVX2 versions would have resulted in a segmentation violation. The SSE2 version
is significantly slower than the AVX512f version, so this situation is best
avoided, though the results will always be correct. Note that the routine
<tt>VectorHashSIMD()</tt> also carries out this test on the buffer alignment, so
if you e.g. instruct <tt>VectorHashSIMD()</tt> to use the AVX2 instruction set,
but supply a 16-byte aligned buffer, the routine will actually use the SSE2
instruction set. On the other hand, <tt>VectorHashSIMD()</tt> will not test the
hardware capabilities of the processor, so this routine can crash if you select
a SIMD instruction set that is not supported by the processor.

In Linux, <tt>malloc()</tt> will typically return buffers aligned on a 16-byte
boundary, so special measures are needed to obtain buffers on a 32- or 64-byte
boundary. The same holds for C++ operators like new or memory allocated in C++
container types as they use the same underlying memory pool. For memory
allocated on the stack even 16-byte alignment is not guaranteed. See [this
page](https://embeddedartistry.com/blog/2017/02/22/generating-aligned-memory/)
for a detailed discussion of how to obtain correctly aligned memory.

### Copyright

VectorHash is distributed with a [zlib open-source software
license](https://opensource.org/licenses/Zlib). You are free to change and
redistribute the software as long as the copyright license is maintained. There
is NO WARRANTY, to the extent permitted by law.

### The cpuid library

VectorHash is distributed with a library for detecting CPU capabilities in the
cpuid subdirectory. This implementation was written by steinwurf. The software
can be found [in this repository](https://github.com/steinwurf/cpuid/). What is
contained in the cpuid subdirectory is a somewhat simplified version of the
9.0.1 release of cpuid. The most important modification is to remove the python
dependency for the build process. The original license and documentation have
been retained, so the latter still refers to the python build process, which no
longer works in this version. The cpuid library is distributed with a BSD-style
copyright license.

### The unittest-cpp library

VectorHash is distributed with a library for unit tests in the unittest-cpp
subdirectory. This implementation was originally written by Noel Llopis and
Charles Nicholson. The software can be found [in this
repository](https://github.com/unittest-cpp/unittest-cpp). What is contained in
the unittest-cpp subdirectory is a somewhat simplified version of the 2.0.0
release. The most important modification is to remove the cmake dependency for
the build process. The original license and documentation have been retained, so
the latter still contains incorrect information about the build process. A
simple Makefile has been added instead. The unittest-cpp library is distributed
with an MIT/X copyright license.


Peter A.M. van Hoof   
Royal Observatory of Belgium   
p.vanhoof@oma.be
