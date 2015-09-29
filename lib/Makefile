#/*--------------------------------------------------------------------------*\
# |  file: Makefile                                                          |
# |                                                                          |
# |  IPOPT MATLAB-interface provided by:                                     |
# |      Enrico Bertolazzi (enrico.bertolazzi@unitn.it)                      |
# |      Dipartimento di Ingegneria Industriale                              |
# |      Universita` degli Studi di Trento                                   |
# |      Via Sommarive 9, I-38123, Trento, Italy                             |
# |                                                                          |
# |  This IPOPT MATLAB-interface is derived from the code by:                |
# |        Peter Carbonetto                                                  |
# |        Dept. of Computer Science                                         |
# |        University of British Columbia, May 19, 2007                      |
# |        Original code is published under the Eclipse Public License.      |
#\*--------------------------------------------------------------------------*/

# INSTRUCTIONS: Please modify the following few variables. See the
# Ipopt documentation (Ipopt/doc/documentation.pdf) for more details.

# This variable corresponds to the base directory of your MATLAB installation.
# This is the directory so that in its 'bin/'  subdirectory you see all the
# matlab executables (such as 'matlab', 'mex', etc.)
#MATLAB_HOME = /Applications/MATLAB_R2015a.app
MATLAB_HOME = /usr/local/MATLAB/R2015b/

# Set the suffix for matlab mex files. The contents of the
# $(MATLAB_HOME)/extern/examples/mex directory might be able to help
# you in your choice.
MEXSUFFIX   = mexmaci64

# This is the command used to call the mex program. Usually, it is
# just "mex". But it also may be something like
# /user/local/R2006b/bin/mex if "mex" doesn't work.
MEX = $(MATLAB_HOME)/bin/mex

#############################################################################
# Do not modify anything below unless you know what you're doing.
#CXX         = clang++
CXX         = g++ -DHAVE_CSTDDEF
CXXFLAGS    = -O3 -DNDEBUG -DIPOPT_BUILD -DMATLAB_MEXFILE

LDFLAGS = $(CXXFLAGS)
INCL    = -I../src  -I/usr/include/coin -I/usr/local/include/coin
LIBS    = -L/usr/local/lib -lipopt -lstdc++

LIBS_STATIC = $(LIBS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W  = echo
MEXFLAGCXX = -cxx -largeArrayDims
MEXFLAGS   = -v $(MEXFLAGCXX) -O CC="$(CXX)" CXX="$(CXX)" LD="$(CXX)"       \
              COPTIMFLAGS="$(CXXFLAGS)" CXXOPTIMFLAGS="$(CXXFLAGS)" \
              LDOPTIMFLAGS="$(LDFLAGS)"

TARGET = ipopt.$(MEXSUFFIX)
SRCS   = ../src/IpoptInterfaceCommon.cc ../src/ipopt.cc

all: $(TARGET)

$(TARGET):
	$(MEX) $(INCL) $(MEXFLAGS) $(LIBS_STATIC) $(SRCS) -output $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCL) -I"$(MATLAB_HOME)/extern/include" -o $@ -c $^

clean:
	rm -f $(OBJS) *.lo $(TARGET)

distclean: clean
