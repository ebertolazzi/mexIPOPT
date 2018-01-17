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

COMPFLAGS  = [ '/O2 ' ] ;
DEFINE     = [ '-DNDEBUG -DMATLAB_MEXFILE -DHAVE_CSTDDEF ' ] ;% -DHAVE_CONFIG_H

%IPOPTBASE = 'C:\Ipopt-3.11.0-Win32-Win64-dll' ;
%IPOPTBASE = 'C:\Ipopt-3.10.1-win32-msvc12_mumps+metis+clapack' ;
%IPOPTBASE = 'C:\Ipopt-3.10.1-Win32-Win64-dll' ;
IPOPTBASE = 'C:\cygwin64\usr\local2' ;


ipoptlib = fullfile(IPOPTBASE,'lib','libipopt.lib');
mumpslib = fullfile(IPOPTBASE,'lib','libcoinmumps.lib');

INCL     = ['-I"' pwd '\..\src" -I"' IPOPTBASE '\include\coin" '] ;

LIBS     = [ '-L' IPOPTBASE '\lib -lipopt -lcoinmumps ' ] ;
MEXFLAGS = [ '-v -cxx -largeArrayDims -outdir .. DEBUGFLAGS="" ' ...
             'COMPFLAGS="$COMPFLAGS' COMPFLAGS '" ' ] ;
cmd = sprintf('mex %s %s %s %s %s',MEXFLAGS,DEFINE,INCL,files,LIBS) ;
disp(cmd);
eval(cmd);
