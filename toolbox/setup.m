old_dir = cd(fileparts(which(mfilename)));

addpath(old_dir);
addpath([old_dir '/lib']);
addpath([old_dir '/bin/windows']);
addpath([old_dir '/bin/windows_mingw']);
addpath([old_dir '/bin/osx']);
addpath([old_dir '/bin/linux_3']);
addpath([old_dir '/bin/linux_4']);
addpath([old_dir '/bin/linux_5']);

cd(old_dir);