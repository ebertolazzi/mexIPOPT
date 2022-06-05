### mexIPOPT with HSL and PARDISO
**by Enrico Bertolazzi**

*IPOPT* can can use HSL library (https://www.hsl.rl.ac.uk) and PARDISO (https://pardiso-project.org).
These libraries can improve the performance of the solver. Due to copyright issues these libraries cannot be included in *IPOPT* executable.
However *IPOPT* can load dynamically the library so that if you have your own copy of these library you can use it with *IPOPT*.

**MA57, MA77, MA86, MA97**

To use these solvers you must have

- `libhsl.dll` on windows
- `libhsl.dylib` on OSX
- `libhsl.so` on LINUX

on a directory accessible by MATLAB or on the same directory of your script that use *IPOPT*.
*HSL* library for *IPOPT* can be obtained from

- https://www.hsl.rl.ac.uk/ipopt/

it is free for academic or personal usage.
To compile the library you may use the CMake script
you find in `cmake_for_ma57_77_86_97` on this repository.
Procedure

- Get the sources file from `https://www.hsl.rl.ac.uk/ipopt/` and unzip
  the files on a dyrectory, e.g. `HSL_FOR_IPOPT`.
- Copy `CMakeLists.txt` and the directory `cmake` into the source file
  directory (e.g. `HSL_FOR_IPOPT`).
- To build the library create directory `build` and inside this directory execute
~~~
cmake ..
make
~~~
- HSL are written in FORTRAN so that for windows OS you must compile the library using MINGW.

If that doesn't work try this:

- https://github.com/coin-or-tools/ThirdParty-HSL

**PARDISO**

You can also use *PARDISO* available at:

- https://pardiso-project.org

pardiso is no more free for personal use, so that you must buy a license to use it.
To use *PARDISO* you must rename (or make an alias)
of the library with

- `libpardiso.dll` on windows
- `libpardiso.dylib` on OSX
- `libpardiso.so` on LINUX

on a directory accessible by MATLAB or on the same directory of your script that use *IPOPT*.

**Author:**

	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
