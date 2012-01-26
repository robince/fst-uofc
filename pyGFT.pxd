#	This software is copyright (c) 2010 UTI Limited Partnership.  
#	The original authors are Robert A. Brown, M. Louis Lauzon 
#	and Richard Frayne.  This software is licensed in the terms 
#	set forth in the "FST License Notice.txt" file, which is 
#	included in the LICENSE directory of this distribution.


cdef extern from "string.h":
	void* memcpy(void *s1, void *s2, size_t n)

cdef extern from "stdlib.h":
	void free(void *ptr)



cdef extern from "gft.h":
	int gft_1dSizeOfPartitions(unsigned int N)
	int* gft_1dPartitions(unsigned int N)
	void gaussian(double *win, int N, int freq)
	void box(double *win, int N, int freq)
	double *windows(int N, void *window)
	void gft_1dComplex64(double *signal, unsigned int N, double *win, int *pars, int stride)
	void gft_2dComplex64(double *image, unsigned int N, unsigned int M, void *window)
	double *gft_1d_interpolateNN(double *signal, unsigned int N, unsigned int M)
