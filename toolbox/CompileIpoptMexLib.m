clc;
clear functions;

old_dir = cd(fileparts(which(mfilename)));

MROOT      = matlabroot;
IPOPT_HOME = '../Ipopt-3.13.2-win64-msvs2019-md';

[~,mexLoaded] = inmem('-completenames');
eval('while mislocked(''ipopt''); munlock(''ipopt''); end;');

disp('---------------------------------------------------------');

CMD = [ 'mex -largeArrayDims -Isrc -I' IPOPT_HOME '/include/coin-or ' IPOPT_HOME '/lib/ipopt.dll.lib ' ] ;
if isunix
  CMD = [CMD, 'CXXFLAGS="\$CXXFLAGS -Wall -O2 -g" '];
elseif ispc
end
if ispc
  CMD = [ CMD, '-output bin/windows/ipopt_win ' ];
elseif ismac
  CMD = [ CMD, '-output bin/osx/ipopt_osx ' ];
elseif isunix
  CMD = [ CMD, '-output bin/linux/ipopt_linux ' ];
else
  error('architecture not supported');
end

CMD = [ CMD, './src/ipopt.cc ./src/IpoptInterfaceCommon.cc ' ];
%    CMD = [ CMD, ...
%      ' CXXFLAGS="\$CXXFLAGS -Wall -O2 -g"' ...
%      ' LDFLAGS="\$LDFLAGS -static-libgcc -static-libstdc++"' ...
%      ' LINKLIBS="-L\$MATLABROOT/bin/\$ARCH -L\$MATLABROOT/extern/bin/\$ARCH -lMatlabDataArray -lmx -lmex -lmat -lm "' ...
%    ];
disp(CMD);
eval(CMD);

cd(old_dir);

disp('----------------------- DONE ----------------------------');
