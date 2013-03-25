/*
 *  gft.h
 *  GFT Framework
 *
 *  Created by Robert Brown on 30/05/08.
 *  Copyright 2008 Robert Brown. All rights reserved.
 *
 */

typedef void (windowFunction)(double*,int,int);

// Windows
void gaussian(double *win, int N, int freq);
void box(double *win, int N, int freq);

// GFT partition and window screen generators
int gft_1dSizeOfPartitions(unsigned int N);
int *gft_1dPartitions(unsigned int N);
int *gft_1dMusicPartitions(unsigned int N, float samplerate, int cents);
double *windows(int N, windowFunction *window);
double *windowsFromPars(int N, windowFunction *window, int *pars);

// 1D GFT Functions
void gft_1dComplex64(double *signal, unsigned int N, double *win, int *pars, int stride);

// 2D GFT Functions
void gft_2dComplex64(double *signal, unsigned int N, unsigned int M, windowFunction *window);

// Interpolation functions
double *gft_1d_interpolateNN(double *signal, unsigned int N, unsigned int M);

// Shift utility
void gft_1d_shift(double *signal, unsigned int N, unsigned int shiftBy);
