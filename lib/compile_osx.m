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

% change this flag to false if mumps do not use MPI support
use_mpi    = true ;
use_static = false ;

% interface sources
files = [ '../src/ipopt.cc ', ...
          '../src/IpoptInterfaceCommon.cc' ] ;

% compiler options
CXXFLAGS = '-fPIC -O3 -DMATLAB_MEXFILE ' ;

% where are headers files
INCL = [ '-I../src ' ...
         '-I/usr/include/coin ' ...
         '-I/usr/local/include ' ...
         '-I/usr/local/include/coin ' ] ;

% libraries for MUMPS
% check if it is installed coinmumps library
if exist('/usr/local/lib/libcoinmumps.a','file') || ...
   exist('/usr/local/lib/libcoinmumps.dylib','file') 
  LIBS = '-L/usr/local/lib -lipopt -lcoinmumps' ;
elseif exist('/usr/local/lib/libmumps_common.dylib','file')
  LIBS = '-L/usr/local/lib -lipopt -lmumps_common -lsmumps -ldmumps -lcmumps -lzmumps -lpord -lscalapack' ;
  % check if non mpi version of mumps is installed
  if exist('/usr/local/lib/libmpiseq.dylib','file')
    use_mpi = false ;
    LIBS    = [ LIBS ' -lmpiseq' ] ;
  end
end

if ~use_mpi
  CXXFLAGS = [ CXXFLAGS ' -DIPOPT_INTERFACE_NO_MPI' ] ;
else
  LIBS = [ LIBS ' -lmpi' ] ;  
end

if use_static
  [status,GFORTRAN] = system('/usr/local/bin/gfortran -print-file-name=libgfortran.a') ;
  [status,QUAD]     = system('/usr/local/bin/gfortran -print-file-name=libquadmath.a') ;
  %[status,GCC]      = system('/usr/local/bin/gfortran -print-file-name=libgcc_ext.10.5.dylib') ;
  [status,GCC]      = system('/usr/local/bin/gfortran -print-file-name=libgcc.a') ;
  files = sprintf('%s %s %s %s ',files, GFORTRAN(1:end-1), QUAD(1:end-1), GCC(1:end-1) ) ;
else
  LIBS = [LIBS ' -L/usr/local/lib/gcc/5 -lgfortran -lquadmath'] ;
end

% frameworks
LIBS2 = '-framework Accelerate ' ;

% compiler options
MEXFLAGS = [ '-v -cxx -largeArrayDims ' ...
             'COPTIMFLAGS="' CXXFLAGS '" ' ...
             'CXXOPTIMFLAGS="' CXXFLAGS '" ' ...
             'LDFLAGS=''$LDFLAGS ' LIBS2 ''''  ] ; % -static -shared-libgcc 

% build and execute compilation command
cmd = sprintf('mex %s %s %s %s',MEXFLAGS,INCL,files,LIBS) ;
eval(cmd);
