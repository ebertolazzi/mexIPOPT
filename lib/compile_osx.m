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
          '../src/IpoptInterfaceCommon.cc ' ] ;

% compiler options
CXXFLAGS = '-fPIC -O3 -DMATLAB_MEXFILE ' ;

% extra libs
LIBS = '' ;

% if you follow the proposed installation instruction 
% the ipopt library is located at /usr/local2
PREFIX = '/opt/local2' ;

INCL = [ '-I../src ' ...
         '-I' PREFIX '/include ' ...
         '-I' PREFIX '/include/coin ' ] ;

% full path where gfortran is installed, MATLAB do not known where is
GFORTRANCMD = '/usr/local/bin/gfortran' ;

% full path MATLAB
MATLAB = ['/Applications/MATLAB_R' version('-release') '.app'] ;

% add list of static version of fortran libraries
[status,GFORTRAN] = system([ GFORTRANCMD ' -print-file-name=libgfortran.a']) ;
[status,QUAD]     = system([ GFORTRANCMD ' -print-file-name=libquadmath.a']) ;
[status,GCC]      = system([ GFORTRANCMD ' -print-file-name=libgcc.a']) ;
files = sprintf('%s %s %s %s',files, ...
                 GFORTRAN(1:end-1), ...
                 QUAD(1:end-1), ...
                 GCC(1:end-1) ) ;

% libraries for MUMPS
% check if it is installed coinmumps library or mumps from homebrew
if exist([ PREFIX '/lib/libcoinmumps.a'],'file') || ...
   exist([ PREFIX '/lib/libcoinmumps.dylib'],'file') 
  LIBS = [ LIBS '-L' PREFIX '/lib -lipopt -lcoinmumps ' ] ;
  % assume that ipopt is compiled with support for ma57 
  % library MA57 inside MATLAB
  LIBS = [ LIBS ' -L' MATLAB '/bin/maci64 -lmwma57 -lmwblas ' ] ;
elseif exist('/usr/local/lib/libmumps_common.dylib','file')
  LIBS = '-L/usr/local/lib -lipopt -lmumps_common -lsmumps -ldmumps -lcmumps -lzmumps -lpord -lscalapack' ;
  % redefine include headers search path
  INCL = '-I../src -I/usr/local/include -I/usr/local/include/coin ' ;

  % check if non mpi version of mumps is installed
  if exist('/usr/local/lib/libmpiseq.dylib','file')
    LIBS     = [ LIBS ' -lmpiseq' ] ;
  else
    LIBS     = [ LIBS ' -lmpi' ] ;
    CXXFLAGS = [ CXXFLAGS ' -DIPOPT_INTERFACE_USE_MPI' ] ;
  end
end

% frameworks
LIBS2 = '-framework Accelerate ' ; % -Wl,-no_compact_unwind

% compiler options -largeArrayDims 
MEXFLAGS = [ '-v -cxx ' ...
             'CXXOPTIMFLAGS="' CXXFLAGS '" ' ...
             'LDFLAGS=''$LDFLAGS ' LIBS2 ''''  ] ;

% build and execute compilation command
cmd = sprintf('mex %s %s %s %s',MEXFLAGS,INCL,files,LIBS) ;
eval(cmd);

disp(cmd);
