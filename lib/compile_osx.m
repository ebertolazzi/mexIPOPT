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

% change this flag if mumps do not use MPI support
use_mpi = true ;

% interface sources
files = [ '../src/ipopt.cc ', ...
          '../src/IpoptInterfaceCommon.cc' ] ;

% compiler options
CXXFLAGS = [ '-fPIC -O3 -DMATLAB_MEXFILE ' ] ;
if use_mpi
  CXXFLAGS = [ CXXFLAGS '-DIPOPT_INTERFACE_USE_MPI ' ] ;
end

% where are headers files
INCL = [ '-I../src ' ...
         '-I/usr/include/coin ' ...
         '-I/usr/local/include ' ...
         '-I/usr/local/include/coin ' ] ;

% libraries linked dynamically
LIBS = [ '-lmumps_common -lsmumps -ldmumps -lcmumps -lzmumps -lpord -lscalapack -lmpi' ] ;

if use_mpi
  LIBS = [ LIBS ' -lscalapack -lmpi' ] ;
else
  LIBS = [ LIBS ' -lmpiseq' ] ;
end

% frameworks and libraries linked statically using
LIBS2 = '/usr/local/lib/libipopt.a -framework Accelerate ' ;

% compiler options
MEXFLAGS = [ '-v -cxx -largeArrayDims ' ...
             'COPTIMFLAGS="' CXXFLAGS '" ' ...
             'CXXOPTIMFLAGS="' CXXFLAGS '" ' ...
             'LDFLAGS=''$LDFLAGS ' LIBS2 ''''  ] ; % -static -shared-libgcc 

% build and execute compilation command
cmd = sprintf('mex %s %s %s %s',MEXFLAGS,INCL,files,LIBS) ;
eval(cmd);
