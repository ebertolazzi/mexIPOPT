### mexIPOPT on LINUX
**by Enrico Bertolazzi**

**How to install IPOPT on LINUX or OSX**

The compilation of a mex file involving external 
shared library is a difficult task cause Matlab
has is internal copy of shared library which
may conflict with the system or user defined
shared library.

In order to produce a working mex file for ipopt
interface the better way is to use static library
as long as you can.

The proposed installation procedure uses /usr/local2 as
installation directory, so that, after generation of
mex file you can delete /usr/local2 ad free memeory.

To make a working static library for IPOPT on my linux
system I di the following procedure:

1) Download the IPOPT tarball from `http://www.coin-or.org/download/source/Ipopt/`

2) Untar/unzip the tarball and go on the root of the distribution. 


notice that metis is disabled. I choose to disable metis and
MPI cause otherwise IPOPT joined with mex interface crash
in special circumstance crash.

7) Configure for compile a static ipopt library:

The following `export` prevent gcc to produce a static
library that cannot be linked with a shared object.
Moreover `FUNNY_MA57_FINT` permits to produce a code
that can link with MA57 solver contained in MATLAB.

~~~
export ADD_CFLAGS="-fPIC -fno-common -DFUNNY_MA57_FINT"
export ADD_CXXFLAGS="-fPIC -fno-common -DFUNNY_MA57_FINT"
export ADD_FFLAGS="-fPIC -fno-common -DFUNNY_MA57_FINT"
~~~

To enable the use of MA57 it is enught modify `config_coinhsl.h.in`
in ThirdParty/HSL/ and change

`#undef COINHSL_HAS_MA57`

to 

`#define COINHSL_HAS_MA57 1`

Undef all the other HSL interfaces unless you have the
corresponding libraries.

The configure script is set to produce **only** a static
library and install the library at `/usr/local2/lib` with
headers at `/usr/local2/include`:

~~~
./configure --enable-static --disable-shared \
            --prefix=/usr/local2 --without-metis \
             -enable-matlab-ma57
make
sudo make install
~~~

after that Ipopt is ready to be linked with the mex interface. 
In `/usr/local2/lib` you shuold find the following library

- libcoinmumps.a
- libipopt.a

To avoid conflict with MATLAB blas/lapack shared
library the mex command must force linking of static version of these libraries. The command is build with the script
`compile_linux.m` or `compile_osx.m`.
Adapt it, if necessary, to meet the configuration of your system.
