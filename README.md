1. Requires MKL, Armadillo and Boost
2. Download Intel MKL and run the .sh file, striaghtforward. Default installation folder is /opt/intel
3. Armadillo 7.8 is compiled and installed at /usr/local
4. Run "source setup_local.sh" to set up LD_LIBRARY_PATH if having error 'cannot open shared object file: No such file or directory'
5. Compile the code as is would generate warnings that the compiler is outdated. Using newer version won't help though, since it is not compatiable with boost and other packages.
