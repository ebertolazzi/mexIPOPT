### mexIPOPT on OSX
**by Enrico Bertolazzi**

**Install IPOPT from source**

To install IPOP from source read `README_how_to_install_ipopt.md`.

**Install IPOPT using homebrew**

To install mexIPOPT on OSX one easy way is to use Homebrew.
Donwload Homebrew from `http://brew.sh` and install.
Then on a shell windows type

~~~
brew install ipopt
~~~

installation may use some options, 
for example to install IPOPT using openblas support:

~~~
brew install ipopt --with-openblas
~~~

for a list of options type

~~~
brew info ipopt
~~~

IPOPT depends on `mumps` which by default use MPI.
It is possible to use `mumps` without MPI support

~~~
brew install mumps --without-mpi
~~~
