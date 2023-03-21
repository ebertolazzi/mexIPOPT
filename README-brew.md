### mexIPOPT and brew
**by Enrico Bertolazzi**

On OSX install ipopt with homebrew `https://brew.sh/index_it` and

~~~
brew install ipopt
brew install openblas
~~~

go into directory bin/osx of the toolbox and run

~~~
ruby copy_dylib.rb
ruby filter_dylib.rb
~~~

in MATLAB run

~~~
>> setup
>> CompileIpoptMexLib
~~~

**Author:**

	Enrico Bertolazzi
	Department of Industrial Engineering
	University of Trento
	enrico.bertolazzi@unitn.it
