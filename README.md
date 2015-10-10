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

The interface is written in C++. 
To compile the mex file

- change the working directory to the `lib` directory.
- run the script `compile_osx` or `compile_linux`

To compile the mex interface you need a static IPOPT library
installed in your system. To install the library read
`README_how_to_install_ipopt.md`.

For the usage of the library use `addpath` to add `lib` directory
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
