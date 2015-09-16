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

files = [ '../src/ipopt.cc ', '../src/IpoptInterfaceCommon.cc' ] ;

CXXFLAGS = [ '-O3 -pipe -DNDEBUG ' ...
             '-Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith ' ...
             '-Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long ' ...
             '-DIPOPT_BUILD -DMATLAB_MEXFILE ' ] ; % -DMWINDEXISINT

LDFLAGS = CXXFLAGS ;

% Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL     = '-I../src -I/usr/local/include/coin'
LIBS     = '-L/usr/local/lib -lipopt -lstdc++'
MEXFLAGS = [ '-v -cxx -largeArrayDims -O ' ...
             'COPTIMFLAGS="' CXXFLAGS '" ' ...
             'CXXOPTIMFLAGS="' CXXFLAGS '" ' ...
             'LDOPTIMFLAGS="' LDFLAGS '"' ] ;

cmd = sprintf('mex %s %s %s %s',MEXFLAGS,INCL,files,LIBS) ;
eval(cmd);
