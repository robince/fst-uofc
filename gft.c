/*
 *  gft.c
 *  GFT Framework
 *
 *  Created by Robert Brown on 30/05/08.
 *	This software is copyright (c) 2010 UTI Limited Partnership.  
 *	The original authors are Robert A. Brown, M. Louis Lauzon 
 *	and Richard Frayne.  This software is licensed in the terms 
 *	set forth in the "FST License Notice.txt" file, which is 
 *	included in the LICENSE directory of this distribution.
 *
 */

#include "gft.h"
#include "fftw3.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.1415926535897931


void fft(int N, double *in, int stride) {
	fftw_plan p;
	p = fftw_plan_many_dft(1,&N, 1, (fftw_complex *)in, NULL, stride, 0, (fftw_complex *)in, NULL, stride, 0, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);	
}


void ifft(int N, double *in, int stride) {
	int i;
	fftw_plan p;
	p = fftw_plan_many_dft(1,&N, 1, (fftw_complex *)in, NULL, stride, 0, (fftw_complex *)in, NULL, stride, 0, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	for (i=0; i<N; i++) {
		in[i*2*stride] /= N+0.0;
		in[i*2*stride+1] /= N+0.0;
	}
}	


void cmul(double *x, double *y) {
	double ac, bd, abcd;
	ac = x[0]*y[0];
	bd = x[1]*y[1];
	abcd = (x[0]+x[1])*(y[0]+y[1]);
	x[0] = ac-bd;
	x[1] = abcd-ac-bd;
}


void cmulByReal(double *x, double multiplier) {
	x[0] *= multiplier;
	x[1] *= multiplier;
}


void shift(double *sig, int N, int amount) {
	double *temp;
	int i,j;
	
	temp = (double *) malloc(N*2*sizeof(double));
	memcpy(temp, sig, N*2*sizeof(double));
	for (i=0; i<N; i++) {
		j = i - amount;
		if (j < 0) j = N-j;
		j = j % N;
		sig[i*2] = temp[j*2];
		sig[i*2+1] = temp[j*2+1];
	}
	free(temp);	
}


void gaussian(double *win, int N, int freq) {
	int i;
	double x;
	double sum;
	for (i = 0; i<N*2; i+=2) {
		x = i / (N*2+0.0);
		win[i] = abs(freq)/sqrt(2*PI)*exp(-pow((x-0.5),2)*pow(freq,2)/2.0);
		win[i+1] = 0.0;
		sum += win[i];
	}
	// Make sure the window area is 1.0
	for (i = 0; i<N*2; i+=2) {
		win[i] /= sum;
	}
	shift(win,N,-N/2);
	fft(N,win,1);
}


void box(double *win, int N, int freq) {
	int i;
	for (i = 0; i < N*2; i+=2) {
		win[i] = 1.0;
		win[i+1] = 0.0;
	}
}


int gft_1dSizeOfPartitions(unsigned int N) {
	return round(log2(N))*2+1;
}


int *gft_1dPartitions(unsigned int N) {
	int sf = 1;
	int ef = 2;
	int cf = 1;
	int width = 1;
	int pcount = 0;
	int pOff;
	int *partitions;
	int sp,ep,sn,en;
	
	partitions = (int *)malloc(sizeof(int)*round(log2(N))*2+1);
	pOff = round(log2(N))*2-1;
	
	while (sf < N/2) {
		sp = cf-width/2-1; ep = cf+width/2-1;
		sn = N-cf-width/2+1; en = N-cf+width/2+1;
		if (ep > N) ep = N;
		if (sn < 0) sn = 0;
		if (width/2 == 0)
			ep+=1; sn-=1;
		partitions[pcount] = ep;
		partitions[pOff-pcount] = en;
		pcount++;
		
		sf = sf+width;
		if (sf > 2) width *= 2;
		ef = sf+width;
		cf = sf+width/2;
	}
	
	partitions[pOff+1] = -1;
	return partitions;
}


int *gft_1dMusicPartitions(unsigned int N, float samplerate, int cents) {
	int i;
	int *partitions;
	float freq;
	float fSpacing = samplerate / (0.0+N);
	float reference = 110.0;
	float logreference = log(reference) / log(2.0);
	float logcent = 1./1200.;
	float logdelta = logcent * cents;
	float max = log(samplerate/2.) / log(2.);
	float min = logreference - (logdelta * floor(logreference / logdelta));
	
	while (pow(2,min+logdelta) - pow(2,min) < fSpacing) {
		min += logdelta;
	}
	
	partitions = (int *) malloc(sizeof(int)*(floor((max-min)/logdelta+2))*2);
	int poff = floor((max-min)/logdelta+2)*2-1;
	for (i = 0; i < (max-min) / logdelta + 1; i++) {
		freq = (min-logdelta/2.) + logdelta*i;
		partitions[i] = round(pow(2,freq) / fSpacing);
		partitions[poff-i] = N-round(pow(2,freq) / fSpacing);
	}

	partitions[poff] = -1;
	return partitions;
}


double *windows(int N, windowFunction *window){
	int fstart, fcentre, fwidth, fend;
	double *fband, *fbandminus;
	double *win, *temp;
	int i;
	
	win = (double *)malloc((N)*sizeof(double)*2);
	temp = (double *)malloc((N)*sizeof(double)*2);
	memset(win, 0, N*sizeof(double)*2);
	
	// For each of the GFT frequency bands.  Using a dyadic scale.

	// Frequency 0 and -1 are special cases
	win[0] = 1.0;
	win[N*2-2] = 1.0;
	
	for (fstart=1; fstart<N/2; fstart*=2) {
		if (fstart < 2) {
			fwidth = 1;
			fcentre = fstart +1;
			fend = fstart+fwidth;
		} else {
			fwidth = fstart;
			fcentre = fstart + fwidth/2 +1;
			fend = fstart + fwidth;
		}
		
		// frequency band that we're working with
		fband = win+fstart*2;
		fbandminus = win+N*2-(fstart*2)-2;
		
		// Construct the window (in the Fourier domain)
		window(temp, N, fcentre);
		shift(temp,N,fwidth/2);
		
		// Put the window in the proper +- partitions
		for (i=0; i<fwidth*2; i+=2) {
			fband[i] = temp[i];
			fband[i+1] = temp[i+1];
			fbandminus[-(i)] = temp[i];
			fbandminus[-(i+1)] = temp[i+1];
		}
	}		
	free(temp);
	return win;
}


double *windowsFromPars(int N, windowFunction *window, int *pars){
	int fstart, fcentre, fwidth, fend;
	double *fband;
	double *win, *temp;
	int i;
	int fcount;
	
	win = (double *)malloc((N)*sizeof(double)*2);
	temp = (double *)malloc((N)*sizeof(double)*2);
	memset(win, 0, N*sizeof(double)*2);
	
	// For each of the GFT frequency bands.  Using a dyadic scale.
	
	// Frequency 0 and -1 are special cases
	win[0] = 1.0;
	win[N*2-2] = 1.0;
	
	
	// For each of the GFT frequency bands
	fcount = 0;
	fstart = 0;
	// frequency band that we're working with
	while (pars[fcount] >= 0) {
		fend = pars[fcount];
		fwidth = fend - fstart;
		if (fstart < N/2) {
			fcentre = fstart + fwidth / 2;
		} else {
			fcentre = abs(-N+fstart) - fwidth / 2;
		}
			
		fband = win+fstart*2;
		if (fend - fstart == 1) {		// Width 1 bands are always weighted 1.0
			fband[0] = 1.0;  fband[0+1] = 0.0;
		} else {
			window(temp, N, fcentre);
			shift(temp,N,fwidth/2);
	
			// Put the window in the proper +- partitions
			for (i=0; i<fwidth*2; i+=2) {
				fband[i] = temp[i];
				fband[i+1] = temp[i+1];
			}
		}			
		fstart = pars[fcount];
		fcount++;
	}
	free(temp);
	return win;
}




void gft_1dComplex64(double *signal, unsigned int N, double *win, int *pars, int stride) {
	int fstart, fend, fcount;
	double *fband;
	int i;
	
	// Do the initial FFT of the signal
	fft(N, signal, stride);
	// Apply the windows	
	for (i=0; i<N*2; i+=2) {
		cmul(signal+i*stride,win+i);
	}
	// For each of the GFT frequency bands
	fcount = 0;
	fstart = 0;
	while (pars[fcount] >= 0) {
		fend = pars[fcount];
		// frequency band that we're working with
		fband = signal+fstart*2*stride;
		// inverse FFT to transform to S-space
		ifft(fend-fstart, fband, stride);
		fstart = pars[fcount];
		fcount++;
	}
}


void gft_2dComplex64(double *image, unsigned int N, unsigned int M, windowFunction *window) {
	int row, col;
	double *win;
	int *pars;
	pars = gft_1dPartitions(N);
	win = windows(N, window);
	for (row=0; row < N; row++) {
		gft_1dComplex64(image+(row*N*2),N,win,pars,1);
	}
	pars = gft_1dPartitions(M);
	win = windows(M, window);
	for (col=0; col < M; col++) {
		gft_1dComplex64(image+(col*2),M,win,pars,N);
	}
	
}


void gft_1d_shift(double *signal, unsigned int N, unsigned int shiftBy) {
	double *temp;
	int i,shiftTo;
	
	temp = (double *)malloc(N*2*sizeof(double));
	memcpy(temp,signal,N*2*sizeof(double));
	for (i = 0; i < N; i++) {
		shiftTo = i+shiftBy;
		if (shiftTo >= N) shiftTo -= N;
		signal[shiftTo*2] = temp[i*2];
		signal[shiftTo*2+1] = temp[i*2+1];
	}
	free(temp);
}


double *gft_1d_interpolateNN(double *signal, unsigned int N, unsigned int M) {
	double *image,*temp;
	int factor;
	int f, t;
	int fstart, fwidth, fcentre, fend, twidth, downBy;
	double *fband, *fbandminus;
	double averageR, averageminusR, averageI, averageminusI;
	int i, iPrime;
	
/*	temp = (double *)malloc(N*2*sizeof(double));
	memcpy(temp,signal,N*2*sizeof(double));
	gft_1d_shift(temp,N,N/2);
	signal = temp; */
	
	if (M == 0) M = N;
	factor = N / (M+0.0);
	image = (double *)malloc(M*M*2*sizeof(double));

	for (fstart=1; fstart<N/2; fstart*=2) {
		if (fstart < 2) {
			fwidth = 1;
			fcentre = fstart +1;
			fend = fstart+fwidth;
		} else {
			fwidth = fstart;
			fcentre = fstart + fwidth/2 +1;
			fend = fstart + fwidth;
		}
		
		// frequency band that we're working with
		fband = signal+fstart*2;
		fbandminus = signal+N*2-(fend*2);

		twidth = M / fwidth;
		downBy = fwidth / M;
			
		for (t=0; t<M; t++) {			
			if (twidth == 0) {						// Downsample
				averageR = averageminusR = 0.0;
				averageI = averageminusI = 0.0;
				for (i=t*downBy-downBy/2; i<t*downBy+downBy/2; i++) {
					if (i < 0) iPrime = fwidth-i-1;
					else iPrime = i;
					averageR += fband[iPrime*2];
					averageI += fband[iPrime*2+1];
					averageminusR += fbandminus[iPrime*2];
					averageminusI += fbandminus[iPrime*2+1];
				}
				image[fstart/factor*M*2+t*2] = averageR;
				image[fstart/factor*M*2+t*2+1] = averageI;
				image[(M-fstart/factor-1)*M*2+t*2] = averageminusR;
				image[(M-fstart/factor-1)*M*2+t*2+1] = averageminusI;
			}
			
			else if (twidth == 1) {					// We're at the right sampling rate already
				image[(fstart/factor)*M*2+t*2] = fband[t*2];
				image[(fstart/factor)*M*2+t*2+1] = fband[t*2+1];
				image[(M-fstart/factor-1)*M*2+t*2] = fbandminus[t*2];
				image[(M-fstart/factor-1)*M*2+t*2+1] = fbandminus[t*2+1];
			}
			
			else {									// Have to interpolate
				i = (t+twidth/2);
				if (i >= M) i = i-M;
				if (i < 0) i = M-i;
				i = i / twidth;
				image[(fstart/factor)*M*2+t*2] = fband[i*2];
				image[(fstart/factor)*M*2+t*2+1] = fband[i*2+1];
				image[(M-fstart/factor-1)*M*2+t*2] = fbandminus[i*2];
				image[(M-fstart/factor-1)*M*2+t*2+1] = fbandminus[i*2+1];
			}
			
			cmulByReal(image+(fstart/factor*M*2+t*2),fwidth);
			cmulByReal(image+((M-fstart/factor-1)*M*2+t*2),fwidth);
			for (f=fstart/factor; f<fend/factor; f++) {
				image[f*M*2+t*2] = image[fstart/factor*M*2+t*2];
				image[f*M*2+t*2+1] = image[fstart/factor*M*2+t*2+1];
				image[(M-f-1)*M*2+t*2] = image[(M-fstart/factor-1)*M*2+t*2];
				image[(M-f-1)*M*2+t*2+1] = image[(M-fstart/factor-1)*M*2+t*2+1];
			}
		}
	}
	/* free(temp); */
	return image;
}

