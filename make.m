% Make script for Matlab

% You need to setup your mexopts to compile
% for GCC I needed to add -std=c99 to allow C++ style comments

% compile GFT c library
mex -c gft.c

% might need to add -I and -L flags to fftw3
cd +gft
mex gft1d.c -I../ ../gft.o -lfftw3
cd ..
