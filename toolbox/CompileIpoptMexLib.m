clc;
clear functions;

old_dir = cd(fileparts(which(mfilename)));

isOctave = isempty(ver('matlab'));

if isOctave % skip for Octave (which doesn't yet have INMEM)
  if ismac
    fprintf('Compile for OCTAVE on MACOSX\n');
  else
    fprintf('Compile for OCTAVE\n');
  end
else
  [~,mexLoaded] = inmem('-completenames');
  eval('while mislocked(''ipopt''); munlock(''ipopt''); end;');
end

disp('---------------------------------------------------------');
SRC = ' ./src/ipopt.cc ./src/IpoptInterfaceCommon.cc ';
CMD = [ 'mex -largeArrayDims -Isrc ' SRC ];

if isOctave
  CMD = [ 'mkoctfile --mex ' SRC ' -Wall -O2 -g -Isrc '];
  if ismac
    %
    % libipopt must be set with:
    % install_name_tool -id "@loader_path/libipopt.3.dylib" libipopt.3.dylib
    %
    IPOPT_HOME = '../Ipopt';
    CMD = [ CMD ...
      '-I' IPOPT_HOME '/include_osx/coin-or ' ...
      '-DOS_OCTAVE_ON_MAC -output bin/osx/ipopt_osx_octave ' ...
      '-Lbin/osx -lcoinmumps.3 -lipopt.3 -lgfortran.5 -lquadmath.0 -lmetis -lstdc++.6 -lgcc_s.1 -ldl ' ...
    ];
  elseif ispc
    IPOPT_HOME = '../Ipopt';
    LIBS = ' "-Wl,-rpath=." -Lbin/windows_mingw ';
    NAMES = {'ipopt','coinmumps'};
    for kkk=1:2
      %LIBS = [ LIBS, ' ./bin/windows_mingw/lib', NAMES{kkk}, '.dll.a ' ];
      LIBS = [ LIBS, ' -l', NAMES{kkk}, '.dll ' ];
    end
    CMD = [ CMD ...
      '-I' IPOPT_HOME '/include_win_mingw/coin-or ' ...
      '-DOS_OCTAVE_ON_WINDOWS -output bin/windows_mingw/ipopt_win_octave ' LIBS ...
    ];
  else
    MEX_EXE = 'bin/ipopt_oct';
    [status, IPOPT_INCL] = system('pkg-config --cflags ipopt', 1);
    [status, IPOPT_LIBS] = system('pkg-config --libs ipopt', 1);
    CMD = [ CMD IPOPT_INCL(1:end-1) ' -output ' MEX_EXE ' ' IPOPT_LIBS(1:end-1) ];
  end
elseif ismac
  %
  % libipopt must be set with:
  % install_name_tool -id "@loader_path/libipopt.3.dylib" libipopt.3.dylib
  %
  %  ' bin/osx/libc++.1.dylib bin/osx/libc++abi.dylib bin/osx/libgcc_s.1.dylib' ...
  IPOPT_HOME = '../Ipopt';
  CMD = [ CMD ...
    '-I' IPOPT_HOME '/include_osx/coin-or ' ...
    '-DOS_MAC -output bin/osx/ipopt_osx ' ...
    ' bin/osx/libcoinmumps.3.dylib bin/osx/libipopt.3.dylib ' ...
    ' bin/osx/libgfortran.5.dylib bin/osx/libquadmath.0.dylib bin/osx/libmetis.dylib ' ...
    ' bin/osx/libstdc++.6.dylib bin/osx/libgcc_s.1.dylib' ...
    'LDFLAGS=''$LDFLAGS -Wl,-rpath,. -framework Accelerate -ldl'' ' ...
    'CXXFLAGS=''$CXXFLAGS -Wall -O2 -g'' ' ...
  ];
  %%  'LDFLAGS=''$LDFLAGS -Wl,-rpath,./ -Lbin/osx -L/usr/local/lib -lipopt -lcoinmumps -lgfortran -lquadmath  -framework Accelerate'' ' ...
  %%  '-output bin/osx/ipopt_osx bin/osx/libgcc_s.1.dylib ' ...
elseif isunix
  IPOPT_HOME = '../Ipopt';
  myCCompiler = mex.getCompilerConfigurations('C','Selected');
  switch myCCompiler.Version
  case {'4','5','6'}
    BIN_DIR = 'bin/linux_3';
    MEX_EXE = 'bin/linux_3/ipopt_linux_3';
  case {'7','8'}
    BIN_DIR = 'bin/linux_4';
    MEX_EXE = 'bin/linux_4/ipopt_linux_4';
  otherwise
    BIN_DIR = 'bin/linux_5';
    MEX_EXE = 'bin/linux_5/ipopt_linux_5';
  end
  CMD = [ CMD ...
    '-I' IPOPT_HOME '/include_linux/coin-or ' ...
    '-DOS_LINUX -output ' MEX_EXE ' '...
    'CXXFLAGS=''$CXXFLAGS -Wall -O2 -g'' ' ...
    'LDFLAGS=''$LDFLAGS -static-libgcc -static-libstdc++'' ' ...
    'LINKLIBS=''-L' BIN_DIR ' -L$MATLABROOT/bin/$ARCH -Wl,-rpath,$MATLABROOT/bin/$ARCH ' ...
              '-Wl,-rpath,. -lipopt -lcoinmumps -lopenblas -lgfortran -lgomp -ldl ' ...
              '-lMatlabDataArray -lmx -lmex -lmat -lm '' ' ...
  ];
elseif ispc
  if true
    % use ipopt precompiled with visual studio
    IPOPT_HOME = '../Ipopt/include_win_mingw/';
    IPOPT_BIN  = 'bin/windows_mingw/';
    LIBS = [' -L' IPOPT_BIN ];
   % NAMES = {'ipopt','sipopt','coinmumps','coinasl','ipoptamplinterface'};
    NAMES = {'ipopt','coinmumps'};
    for kkk=1:length(NAMES)
      LIBS = [ LIBS, ' -l', NAMES{kkk} ];
    end
    CMD = [ CMD ...
      '-DOS_WIN -I' IPOPT_HOME '/coin-or ' ...
      '-output ' IPOPT_BIN 'ipopt_win_mingw ' LIBS ...
    ];
  else
    % use ipopt precompiled with visual studio
    IPOPT_HOME = '../Ipopt/include_vs/';
    IPOPT_BIN  = 'bin/windows/';
    LIBS = [' -L' IPOPT_BIN ];
    NAMES = {'ipoptfort','ipopt','libomp','ompstub','openblas','flang','flangmain','flangrti'};
    for kkk=1:8
      LIBS = [ LIBS, ' -l', NAMES{kkk} ];
    end
    CMD = [ CMD ...
      '-DOS_WIN -I' IPOPT_HOME '/coin-or ' ...
      '-output ' IPOPT_BIN 'ipopt_win ' LIBS ...
    ];
  end
else
  error('architecture not supported');
end

disp(CMD);
eval(CMD);

cd(old_dir);

disp('----------------------- DONE ----------------------------');
