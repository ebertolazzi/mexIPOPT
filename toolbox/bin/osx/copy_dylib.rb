require 'colorize'
#
# Adjust path for depend on dylib in the local bin directory
#
require 'fileutils'

fortran_home = Dir["/usr/local/Cellar/gcc/*/lib/gcc/current"].last
ipopt_hdr    = Dir["/usr/local/Cellar/ipopt/*/include"].last
ipopt_lib    = Dir["/usr/local/Cellar/ipopt/*/lib"].last
openblas_lib = Dir["/usr/local/Cellar/openblas/*/lib"].last

FileUtils.rm_rf "../../src_ipopt_osx"
FileUtils.mkdir_p "../../src_ipopt_osx"
FileUtils.cp_r "#{ipopt_hdr}/coin-or", "../../src_ipopt_osx"

[
  'libgfortran.5.dylib',
  'libquadmath.0.dylib',
  'libstdc++.6.dylib',
  'libgcc_s.1.dylib'
].each do |lib|
  #puts fortran_home+"/"+lib, lib
  FileUtils.rm lib if File.exist? lib
  FileUtils.cp fortran_home+"/"+lib, lib
end

[
  'libopenblas.0.dylib',
  'libblas.dylib',
  'liblapack.dylib'
].each do |lib|
  #puts fortran_home+"/"+lib, lib
  FileUtils.rm lib if File.exist? lib
  FileUtils.cp "#{openblas_lib}/#{lib}", lib
end

[
  'libdmumps.dylib',
  'libipopt.3.dylib',
  'libipoptamplinterface.3.dylib',
  'libmpiseq.dylib',
  'libmumps_common.dylib',
  'libpord.dylib',
  'libsipopt.3.dylib'
].each do |lib|
  FileUtils.rm lib if File.exist? lib
  FileUtils.cp ipopt_lib+"/"+lib, lib
end
