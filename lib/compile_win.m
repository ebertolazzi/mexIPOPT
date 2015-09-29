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

CXXFLAGS = [ '/O3 /DNDEBUG /DIPOPT_BUILD /DMATLAB_MEXFILE /DHAVE_CSTDDEF' ] ; % -DMWINDEXISINT

LDFLAGS = CXXFLAGS ;

%IPOPTBASE = 'C:/Ipopt-3.10.1-win32-msvc10_mumps+metis+clapack'
IPOPTBASE = 'C:/Ipopt-3.11.1-win32-cl16ifort13.1'
%IPOPTBASE = '"C:/Program Files (x86)/COIN-OR/1.7.4/win32-msvc10"'

INCL     = ['-I../src  -I' IPOPTBASE '/include/coin '] ;
LIBS     = ['-L' IPOPTBASE '/lib -lipopt ' ] ;
MEXFLAGS = [ '-v -cxx -largeArrayDims -outdir .. ' ...
             'COPTIMFLAGS="' CXXFLAGS '" ' ...
             'CXXOPTIMFLAGS="' CXXFLAGS '" ' ...
             'LDOPTIMFLAGS="' LDFLAGS '"' ] ;
cmd = sprintf('mex %s %s %s %s',MEXFLAGS,INCL,files,LIBS) ;
eval(cmd);
