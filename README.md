### mexIPOPT
**by Enrico Bertolazzi**

This is my rewrite of Peter Carbonetto MATLAB interface for IPOPT
a software package for large-scale ​nonlinear optimization.

Source code and documentation for IPOPT can be downloaded from:

`https://projects.coin-or.org/Ipopt`

On OSX IPOPT can be installed using Homebrew (`http://brew.sh`).

**Why do a job already done?**

The original MATLAB interface (`https://projects.coin-or.org/Ipopt/wiki/MatlabInterface`) seems not maintained for OSX and do not work with recent MATLAB distributions.
I have reorganized and simplified the original interface (only internally) and eliminated the bug (due to MATLAB change 
in the managing of sparse pattern) in the interfacing
of sparse pattern between MATLAB and IPOPT.
Moreover I tried to improve error catching.

**How to install**

Use the toobox installer at

- `https://github.com/ebertolazzi/mexIPOPT/releases`

download and run in MATLAB to install following the instruction.

**How to compile**

Compilation shuold be not necessary. To (re-)compile the mex file
for your architecture

- change the working directory to the `toolbox` directory.
- run the script `CompileIpoptMexLib`.

if all the thing go well open 

- run the script `setup` and them `../IPOPT-toolbox.prj` to compile the toolbox.

at this point you find the recompiled toolbox at

- `../IPOPT-toolbox.mltbx`

**Examples**

In the directory `toolbox/example` you find the original examples 
of Peter Carbonetto which shows the usage of the interface.

In the directory `test_TRAIN` you find a complex example
of an Optimal Control Problem solved using direct trascription
and transformed to an NLP solved using IPOPT.
The descriptin of the problem is at

- `https://vanderbei.princeton.edu/tex/trajopt/trajopt.pdf`.

**Author:**
	
	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
