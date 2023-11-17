### mexIPOPT and brew
**by Enrico Bertolazzi**

On OSX install ipopt with homebrew `https://brew.sh/index_it` and

~~~
brew install ipopt
brew install openblas
brew install metis
~~~

after installation compile MEX with cmake

~~~
cd toolbox
cmake -Bbuild .
cmake --build build --config Release
cmake --install build --config Release
~~~

if all the things are done well you find `ipopt.mexmaci64`
in the directory `lib` ready to be used.


Portable binary
---------------

The compiled mex depend on the presence in the system of ipopt installed with homebrew.
To avoid this dependence and produce a portable toolbox copy the file `copy_dylib.rb` and `filter_dylib.rb` from
the direftory `precompiled_binary/osx` into `toolbox/lib` where is istalled `ipopt.mexmaci64`. Run both script

~~~
ruby copy_dylib.rb
ruby filter_dylib.rb
~~~

if all the things go well now the directory is populated with
the dependent library and `rpath` of the libraries are
changes reflecting the new position.

**Author:**

	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
