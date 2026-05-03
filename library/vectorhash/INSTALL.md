## VectorHash - a very fast vectorized hash function

### Compiling the code

A configuration step is not needed, so we can go straight to the actual
compilation. A GNU Makefile is supplied, so you can do the compilation simply by
typing:

	make

On some operating systems you may have to use <tt>gmake</tt> instead of
<tt>make</tt>. If you want to do a parallel build, e.g. using 16 cores, you can
type:

	make -j16

The build will create three executables called vh128sum, vh256sum, and
vh512sum. These are hard-linked to each other, so they really only take up the
space of one executable. The Makefile is set up to use the GNU compiler g++.
It is straightforward to use the LLVM compiler clang++ instead:

	make -j16 CXX=clang++

In principle it should be possible to use other C++11 compatible compilers too,
but this has not been tested. In that case you will most likely have to change
the compiler options as well. Assuming your favorite compiler is called
<tt>CC</tt>, you would get something like this:

	make -j16 CXX=CC CXXFLAGS="CC-specific-flags"

The build process will also create a 64-bit shared object library called
<tt>libvhsum.so</tt>. This can be used to checksum buffers in your own programs.
See the file <tt>README.md</tt> for more details.

#### 32-bit systems

The text above assumed that you are compiling the code on a 64-bit operating
system. These are now ubiquitous. However, limited support for legacy 32-bit
systems is available. I do not have access to a native 32-bit system, but it
is expected that the build as described above will work on such systems as
well. You will obviously get 32-bit executables and libraries instead of the
64-bit versions mentioned above. On such systems the library will internally
be stored in the lib64 directory, but will eventually be installed in the
correct directory.

#### 32-bit library

On 64-bit systems, the default build will only create a 64-bit shared object
library. If you would like to have an additional 32-bit shared object library,
you create this with the command (after the previous step is complete):

	make -j16 lib32

### Testing the code

After the build is complete, you can optionally test the binaries and the
library by issuing the command:

	make -j16 check

If you created a 32-bit library on a 64-bit operating system, you can test
that library as well by using these commands:

	make testclean
	make -j16 check32

Note that this second step is not needed if you built the library on a 32-bit
operating system. The first step will already have tested that library.

### Installing the code

Installing the code can be done by simply typing:

	make install

You may need root privileges to install the files in the location you have
chosen. In that case you will need to log in as root or use the <tt>sudo</tt>
command:

	sudo make install

By default the installation will be done in <tt>/usr/local</tt>. If you want to
use a different directory, you can type:

	make install INSTALLDIR=/path/to/dir

This will install files in <tt>INSTALLDIR/bin</tt> (for the executables),
<tt>INSTALLDIR/lib*</tt> (for the libraries, more on that below),
<tt>INSTALLDIR/include</tt> (for the header file), and
<tt>INSTALLDIR/man/man*</tt> (for the man pages). The default location for
installing 64-bit libraries varies for different operating systems. For instance
on macOS or Debian-based Linux distros, they would need to be installed in
<tt>/usr/local/lib</tt> while on other Linux distros they would need to be in
<tt>/usr/local/lib64</tt>. The build system has heuristics built in that will
determine this automatically.

Peter A.M. van Hoof   
Royal Observatory of Belgium   
p.vanhoof@oma.be
