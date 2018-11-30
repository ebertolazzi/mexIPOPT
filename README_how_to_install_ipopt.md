### IPOPT compilation
**by Enrico Bertolazzi**

**How to compile IPOPT on LINUX or OSX**

The compilation of a mex file involving external shared library is a
difficult task cause Matlab has is internal copy of shared library 
which may conflict with the system or user defined
shared library.

In order to produce a working mex file for ipopt interface the better
way is to use static library as long as you can.

The proposed installation procedure uses `/opt/local2`
as installation directory, so that, after generation of
mex file you can delete /usr/local2 and free the memory.

To make a working static library for IPOPT on my linux
system I done the following procedure:

1) Download the IPOPT tarball from `http://www.coin-or.org/download/source/Ipopt/`

2) Untar/unzip the tarball and go on the root of the distribution. 

3) get Mumps with the commands:

~~~
cd ThirdParty/Mumps
./get.Mumps
~~~

4) Configure for compile a static mumps library:

~~~
export ADD_CFLAGS="-fPIC -fno-common -mmacosx-version-min=10.9"
export ADD_CXXFLAGS="-fPIC -fno-common -mmacosx-version-min=10.9"
export ADD_FFLAGS="-fPIC -fno-common -mmacosx-version-min=10.9"
./configure --prefix=/opt/local2 --enable-static --disable-shared --without-metis
make
make install
~~~

(note directory /opt/local2 can be written with the new OS)

7) Configure for compile a static ipopt library:

The following `export` prevent gcc to produce a static
library that cannot be linked with a shared object.
Moreover `FUNNY_MA57_FINT` permits to produce a code
that can link with MA57 solver contained in MATLAB.

~~~
export ADD_CFLAGS="-fPIC -fno-common -DFUNNY_MA57_FINT -DHAVE_STDIO -DHAVE_STRING -mmacosx-version-min=10.9"
export ADD_CXXFLAGS="-fPIC -fno-common -DFUNNY_MA57_FINT -DHAVE_STDIO -DHAVE_STRING -mmacosx-version-min=10.9"
export ADD_FFLAGS="-fPIC -fno-common -DFUNNY_MA57_FINT -DHAVE_STDIO -DHAVE_STRING -mmacosx-version-min=10.9"
~~~

To enable the use of MA57 it is enught modify `config_coinhsl.h.in`
in ThirdParty/HSL/ and change

`#undef COINHSL_HAS_MA57`

to 

`#define COINHSL_HAS_MA57 1`

Undef all the other HSL interfaces unless you have the
corresponding libraries.
Change also in `ipopt/src/contrib/LinearSolverLoader/HSLLoader.cc` the line with

~~~
#define HSLLIBNAME "libhsl." SHAREDLIBEXT
~~~

with

~~~
#define HSLLIBNAME "libmwma57." SHAREDLIBEXT
~~~

otherwise when load `ma57` the mex try to load the dynamic library `libhsl.EXT` where EXT = {dylib,so,dll} instead of the MATLAB version.

The configure script must be set to produce **only** a static
library and install the library at `/opt/local2/lib` with
headers at `/opt/local2/include`:

~~~
./configure --enable-static --disable-shared --prefix=/opt/local2 --without-metis --enable-matlab-ma57
~~~

notice that metis is disabled. I choose to disable metis and
MPI otherwise IPOPT joined with mex interface crash
in special circumstance.

**On LINUX**
On linux OSX if you want to use MATLAB `ma57` you must compile Blas/Lapack and use its static version

~~~
cd ThirdParty/Blas
./get.Blas
cd ../../ThirdParty/Lapack
./get.Lapack
~~~

moreover the configure scripts must be changed as 

~~~
export ADD_CFLAGS="-fpic -fno-common -DFUNNY_MA57_FINT -DHAVE_STDIO -DHAVE_STRING"
export ADD_CXXFLAGS="-fpic -fno-common -DFUNNY_MA57_FINT -DHAVE_STDIO -DHAVE_STRING"
export ADD_FFLAGS="-fpic -fno-common -DFUNNY_MA57_FINT -DHAVE_STDIO -DHAVE_STRING"
~~~

and

~~~
./configure --enable-static --disable-shared --prefix=/opt/local2 --without-metis --enable-matlab-ma57 --with-pic --with-blas=BUILD --with-lapack=BUILD 
~~~

**Compile IPOPT**

~~~
make
sudo make install
~~~

after that Ipopt is ready to be linked with the mex interface. 
In `/usr/local2/lib` you should find the following library

- libcoinmumps.a
- libipopt.a

To avoid conflict with MATLAB blas/lapack shared
library the mex command must force linking of 
static version of these libraries.
The command is build with the script

`compile_linux.m` or `compile_osx.m`.

Adapt it, if necessary, to meet the configuration of your system.

**PARDISO NOTE**

Also pardiso solver ('http://www.pardiso-project.org') can be used in alternative to `ma57`.
Pardiso can be activated but normally libraries are linked with dylib located at `/usr/local/lib/` while with the new MacOS libraries are located at different location, for example at `/usr/local/Cellar/gcc/7.3.0/lib/gcc/7/`. In this case change the library search path using the command `install_name_tool` for example:

~~~
install_name_tool -change /usr/local/lib/libgfortran.3.dylib /usr/local/Cellar/gcc/7.3.0/lib/gcc/7/libgfortran.dylib libpardiso.dylib

install_name_tool -change /usr/local/lib/libgomp.1.dylib /usr/local/Cellar/gcc/7.3.0/lib/gcc/7/libgomp.dylib libpardiso.dylib

install_name_tool -change /usr/local/lib/libquadmath.0.dylib /usr/local/Cellar/gcc/7.3.0/lib/gcc/7/libquadmath.dylib libpardiso.dylib
~~~
