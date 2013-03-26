% Make script for Matlab

% You need to setup your mexopts to compile
% for GCC I needed to add -std=c99 to allow C++ style comments

% compile GFT c library
% here need -I to point to fftw headers
mex -c gft.c -I/usr/local/include

% here need -L to point to fftw3 library
cd +gft
mex gft1d.c -I../ ../gft.o -lfftw3 -L/usr/local/lib
mex gft1dInterpolate.c -I../ ../gft.o -lfftw3 -L/usr/local/lib
mex gft2d.c -I../ ../gft.o -lfftw3 -L/usr/local/lib
cd ..
