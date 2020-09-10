### mexIPOPT
**by Enrico Bertolazzi**

This is my rewrite of Peter Carbonetto MATLAB interface for IPOPT
a software package for large-scale ​nonlinear optimization.

Source code and documentation for IPOPT can be downloaded from:

`https://projects.coin-or.org/Ipopt`

On OSX IPOPT can be installed using Homebrew (`http://brew.sh`).

**Why do a job already done?**

The original MATLAB interface (`https://projects.coin-or.org/Ipopt/wiki/MatlabInterface`) seems not maintained for OSX and do not 
work with recent MATLAB distributions.
I have reorganized and simplified the original interface 
(only internally) and eliminated the bug (due to MATLAB change 
in the managing of sparse pattern) in the interfacing
of sparse pattern between MATLAB and IPOPT.
Moreover I tried to improve error catching.

**How to compile**

If you have an OSX or a LINUX 64 bit OS try to use
the corresponding precompiled mex file, otherwise
try to compile the mex interface following the 
instructions.
Compilation is not an easy task, I try to be
as clear as possibile, in any case things may
change from different version of OS.

The interface is written in C++.
To compile the mex file

- change the working directory to the `lib` directory.
- run the script `compile_osx` or `compile_linux` or `compile_win`. 

To compile the mex interface you need a valid IPOPT library
installed in your system. To install the library read
`README_osx_how_to_compile_ipopt.md` or
`README_linux_how_to_compile_ipopt.md`.
For windows system in the `binary` directory I copied precompiled
dll of ipopt + dependecties taken from `https://github.com/JuliaOpt`.

For the usage of the mex file use `addpath` to add `lib` directory
in the search path or move the contents of `lib` to a directory 
in the search path.

**Examples**

In the directory `example` you find the original examples 
of Peter Carbonetto which shows the usage of the interface.

**Author:**
	
	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
