### mexIPOPT
**by Enrico Bertolazzi**

This is my rewrite of Peter Carbonetto MATLAB interface for IPOPT
a software package for large-scale ​nonlinear optimization.

Source code and documentation for IPOPT can be downloaded from:

`https://projects.coin-or.org/Ipopt`

On OSX IPOPT can be installed using Homebrew (`http://brew.sh`).

**Why do a job already done?**

The original MATLAB interface (`https://projects.coin-or.org/Ipopt/wiki/MatlabInterface`) seems not maintained for OSX and do not work with recent
MATLAB distributions.
I have reorganized and simplified the original interface 
(only internally) and eliminated the bug in the interfacing
of sparse pattern between MATLAB and IPOPT.
Moreover I tried to improve error catching.

**How to compile**

The interface is written in C++. 
To compile the mex file

- change the working directory to the `lib` directory.
- run the script `compile_osx` or `compile_linux`

For the uage of the library use `addpath` to add `lib` directory
in the search path or move the contents of `lib` to a directory 
in the search path.

**Note on static library**

MATLAB uses it own copy of dynamic libraries. In particular
there may be conflict with `libgfortran.3.dylib`.
If it happen try to link with the static version of
the fortran libraries by setting:

~~~
use_static = true ;
~~~

in `compile_osx` blas/lapack libraries are included 
in the system framework `Accelerate` and are not
statically linked.
Linux version are statically linked with
fortran and blas/lapack libraries.

**Note on MPI**

The IPOPT library normally interfaces with `MUMPS` which need
initialization of Message Passing Interface (MPI).
This also imply that the interface need to be interfaced with 
MPI library.
If IPOPT is not compiled using MPI library you may remove
the MPI initialization on the mex interface by changing

~~~
use_mpi = false ;
~~~

on the script `compile_osx.m` or `compile_linux.m`
which define the macro `IPOPT_INTERFACE_NO_MPI` and 
remove MPI support on the interface.

**Examples**

In the directory `example` you find the original examples 
of Peter Carbonetto which shows the usage of the interface.

**Author:**
	
	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
