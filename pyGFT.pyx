#	This software is copyright (c) 2010 UTI Limited Partnership.  
#	The original authors are Robert A. Brown, M. Louis Lauzon 
#	and Richard Frayne.  This software is licensed in the terms 
#	set forth in the "FST License Notice.txt" file, which is 
#	included in the LICENSE directory of this distribution.


from pyGFT cimport *
import numpy as np
cimport numpy as np

ntype = np.complex128
ctypedef np.complex128_t ntype_t


def phase(num):
	try:
		num = num.astype(ntype)
	except:
		pass
	angle = np.arctan2(num.imag,num.real)
	return angle
	
	
def complexWithRP(r,p):
	return np.complex(r*np.cos(p),r*np.sin(p))

def shift2D(a,nx,ny):
	b = np.concatenate((a[ny:,:],a[0:ny,:]))
	c = np.concatenate((b[:,nx:],b[:,0:nx]),1)
	return c

def shift1D(a,nx):
	nx = -int(round(nx))
	b = np.concatenate((a[nx:],a[0:nx]))
	return b
	

def shift(a,nx,ny=None):
	if len(np.shape(a)) == 1:
		b = shift1D(a,nx)
	else:
		b = shift2D(a,nx,ny)
	return b


def partitions(N):
	cdef int *a
	cdef int n = gft_1dSizeOfPartitions(N)
	cdef np.ndarray arr
	try:
		if (N.__class__ == int().__class__):
			a = gft_1dPartitions(N)
		else:
			a = gft_1dPartitions((np.max(np.shape(N))))
	except:
		a = gft_1dPartitions((np.max(np.shape(N))))
	arr = np.zeros(n,dtype=np.int)
	memcpy(<np.int_t*>arr.data,a,sizeof(int)*n)
	arr = arr.compress(np.where(arr>0,1,0))
	free(a)
	return arr



def gaussianWindows(N):
	cdef np.ndarray sig
	cdef void *window = &gaussian
	cdef double *win
	win = windows(N,window)
	sig = np.zeros(N,ntype)
	memcpy(<ntype_t*>sig.data,win,N*2*sizeof(double))
	free(win)
	return sig.real


def gft1d(sig,windowType='gaussian'):
	sig = sig.astype(ntype)
	cdef int N = len(sig)
	cdef void *window
	cdef int *pars = gft_1dPartitions(N)
	cdef double *win
	cdef np.ndarray ret = sig.copy()
	if (windowType == 'gaussian'):
		window = &gaussian
	else:
		window = &box
	win = windows(N,window)
	gft_1dComplex64(<double *>ret.data, N, win, pars, 1);		
	ret = shift(ret,N/2)
	free(pars)
	free(win)
	return ret
	

def gft2d(image,windowType='gaussian'):
	cdef int N, M
	(N,M) = np.shape(image)
	image = image.astype(ntype)	
	cdef void *window
	cdef np.ndarray ret = image.copy()
	if (windowType == 'gaussian'):
		window = &gaussian
	else:
		window = &box
	gft_2dComplex64(<double *>ret.data, N, M, window);		
	ret = shift2D(ret,N/2,M/2)
	return ret


def gft1dInterpolateNN(SIG,M=None):
	cdef double *image
	cdef np.ndarray ret
	cdef np.ndarray s
	N = len(SIG)
	if not M:
		M = N
	ret = np.zeros((M,M),ntype)
	s = shift(SIG.astype(ntype),-N/2)
	image = gft_1d_interpolateNN(<double*>s.data,N,M)
	memcpy(<ntype_t*>ret.data,image,M*M*2*sizeof(double))
	ret = np.reshape(ret,(M,M))
	free(image)
	return ret
