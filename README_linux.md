### mexIPOPT on LINUX
**by Enrico Bertolazzi**

**How to install IPOPT**

The compilation of a mex file involving external shared
library is difficult task cause Matlab has is inetrnal 
copy of shared library which often conflict with the system
or user defined shared library.
In order to produce a working mex file for ipopt interface
the better way is to use static library as long as you can.

To make a working static library for IPOPT on my linux
system I di the following procedure:

1) Download the IPOPT tarball from `http://www.coin-or.org/download/source/Ipopt/`

2) Untar/unzip the tarball and go on the root of the distribution. Then get mumps with the commands:

~~~
cd ThirdParty/Mumps
./get.Mumps
~~~

3) Configure for compile a static ipopt linbrary:

The following `export` prevent gcc to produce a static
library that cannot be linked with a shared object:

~~~
export ADD_CFLAGS="-fPIC -fno-common"
export ADD_CXXFLAGS="-fPIC -fno-common"
export ADD_FFLAGS="-fPIC -fno-common"
~~~

The configure script is set to produce **only** a static
library and install the library at `/usr/local/lib` with
headers at `/usr/local/include`:

~~~
./configure --prefix=/usr/local --enable-static --disable-shared

make
make install
~~~

after that the library is ready to be linked with the Ipopt
interface. To avoid conflict with MATLAB blas/lapack shared
library the mex command must force linking of static version of these libraries. The command is build with the script
`compile_linux.m`. Adapt it if necessary to meet the configuration of your system.


