%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  file: Compile.m                                                         %
%                                                                          %
%  IPOPT MATLAB-interface provided by:                                     %
%      Enrico Bertolazzi (enrico.bertolazzi@unitn.it)                      %
%      Dipartimento di Ingegneria Industriale                              %
%      Universita` degli Studi di Trento                                   %
%      Via Sommarive 9, I-38123, Trento, Italy                             %
%                                                                          %
%  This IPOPT MATLAB-interface is derived from the code by:                %
%        Peter Carbonetto                                                  %
%        Dept. of Computer Science                                         %
%        University of British Columbia, May 19, 2007                      %
%        Original code is published under the Eclipse Public License.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interface sources
files = [ '../src/ipopt.cc ', ...
          '../src/IpoptInterfaceCommon.cc' ] ;

% compiler options
CXXFLAGS = [ '-fPIC -O3 -DMATLAB_MEXFILE -DHAVE_CSTDDEF' ] ; % -DMWINDEXISINT

% where are headers files
INCL = [ '-I../src ' ...
         '-I/usr/include/coin ' ...
         '-I/usr/local/include/coin ' ] ;

% libraries linked dynamically
LIBS     = [ '-L/usr/local/lib -L/usr/lib -L/usr/lib/openblas-base ' ...
              '-lgfortran -lquadmath -ldl -lstdc++ -lgcc -lgcc_s -lm -lc' ] ;

% libraries linked statically using
LIBS2    = '-Wl,-Bstatic -lipopt -lcoinmumps -llapack -lopenblas -Wl,-Bdynamic ';

% compiler options
MEXFLAGS = [ '-v -cxx -largeArrayDims ' ...
             'COPTIMFLAGS="' CXXFLAGS '" ' ...
             'CXXOPTIMFLAGS="' CXXFLAGS '" ' ...
             'LDOPTIMFLAGS="" ' ...
             'CXXLIBS=''$CXXLIBS ' LIBS2 ''''  ] ; % -static -shared-libgcc 

% build and execute compilation command
cmd = sprintf('mex %s %s %s %s',INCL,files,MEXFLAGS,LIBS) ;
eval(cmd);
