old_dir = cd(fileparts(which(mfilename)));

addpath('.');
addpath('./lib');
addpath('./bin/windows');
addpath('./bin/osx');
addpath('./bin/linux_3');
addpath('./bin/linux_4');
addpath('./bin/linux_5');

cd(old_dir);