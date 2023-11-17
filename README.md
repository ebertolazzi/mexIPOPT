[![View ebertolazzi/mexIPOPT on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/53040-ebertolazzi-mexipopt)

### mexIPOPT
**by Enrico Bertolazzi**

This is my rewrite of Peter Carbonetto MATLAB interface for IPOPT
a software package for large-scale ​nonlinear optimization.

Source code and documentation for IPOPT can be downloaded from:

https://projects.coin-or.org/Ipopt

On OSX IPOPT can be installed using Homebrew (`http://brew.sh`).

**Why do a job already done?**

The original MATLAB interface (https://projects.coin-or.org/Ipopt/wiki/MatlabInterface) seems not maintained for OSX and do not work with recent MATLAB distributions.
I have reorganized and simplified the original interface (only internally) and eliminated the bug (due to MATLAB change
in the managing of sparse pattern) in the interfacing
of sparse pattern between MATLAB and IPOPT.
Moreover I tried to improve error catching.

**How to install**

Use the toobox installer at

- https://github.com/ebertolazzi/mexIPOPT/releases

**Examples**

In the directory `toolbox/example` you find the original examples
of Peter Carbonetto which shows the usage of the interface.

In the directory `test_TRAIN` you find a complex example
of an Optimal Control Problem solved using direct trascription
and transformed to an NLP solved using IPOPT.
The descriptin of the problem is at

- https://vanderbei.princeton.edu/tex/trajopt/trajopt.pdf.

**MA57, MA77, MA86, MA97, PARDISO**

*IPOPT* can can use HSL library (https://www.hsl.rl.ac.uk) and PARDISO (https://pardiso-project.org).
These libraries can improve the performance of the solver. Due to copyright issues these libraries cannot be included in *IPOPT* executable. Read README_HSL.md to see how to use these libraries in your applications.

**Author:**

	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
