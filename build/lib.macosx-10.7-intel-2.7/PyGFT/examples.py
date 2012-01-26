""" Examples of the 1D and 2D GFT.  Take a look at simpleExample()."""

#	This software is copyright (c) 2010 UTI Limited Partnership.
#	The original authors are Robert A. Brown, M. Louis Lauzon
#	and Richard Frayne.  This software is licensed in the terms
#	set forth in the "FST License Notice.txt" file, which is
#	included in the LICENSE directory of this distribution.

# import pylab for plotting, and PyGFT for doing the GFT
import pylab
import numpy
import PyGFT as gft

# Set up pylab defaults to make the plots prettier
pylab.rcParams['image.origin'] = 'lower'
pylab.rcParams['figure.subplot.bottom'] = 0.07
pylab.rcParams['figure.subplot.top'] = 0.95
pylab.rcParams['figure.subplot.right'] = 0.97
pylab.rcParams['figure.subplot.left'] = 0.06
pylab.rcParams['figure.subplot.wspace'] = 0.10
pylab.rcParams['figure.subplot.hspace'] = 0.10



# Some utility functions


# Produce an image where the pixel values are the distance from one corner
def distCorner(a):
	s = numpy.shape(a)
	for x in range(0,s[0]):
		for y in range(0,s[1]):
			a[x,y] = numpy.sqrt(x**2+y**2)
	return a

# Produce an image where the pixel values are the distance from the centre
def dist(a):
	s = numpy.shape(a)
	cx = s[0] / 2.
	cy = s[0] / 2.
	for x in range(0,s[0]):
		for y in range(0,s[1]):
			a[x,y] = numpy.sqrt((cx-x)**2+(cy-y)**2)
	return a


# Produce a Gaussian
def gaussianWindow(signal,t,k):
	l = len(signal)
	x = numpy.arange(0,1,1./l).astype(numpy.complex64)
	win = 1*numpy.abs(k)/numpy.sqrt(2*numpy.pi)*numpy.e**(-(x-0.5)**2*numpy.abs(k)**2/2.)
	win = win/(numpy.sum(win)/l)
	win = gft.shift(win,-l/2)
	win = gft.shift(win,t)
	return win


# Produce a 2D, circularly symmetric Gaussian
def gaussian2D(image,k):
	shp = numpy.shape(image)
	N = shp[0]
	win = gaussianWindow(image[0],N/2,k)
	win = numpy.array(numpy.mat(win).transpose() * numpy.mat(win))
	return win


# Plot a 1D signal and it's GFT, with the windows and partitions
def plot1D(signal,window,figTitle=None):
	N = len(signal)
	x = numpy.arange(0,1,1./N)
	transform = gft.gft1d(signal,window)
	pylab.figure()
	pylab.subplot(211)
	if figTitle: pylab.title(figTitle)
	pylab.plot(x,signal.real,'b')
	pylab.plot(x,signal.imag,'b--')
	ax = pylab.axis()

	pylab.subplot(212)
	pylab.plot(x,transform.real,'b')
	pylab.plot(x,transform.imag,'g--')
	partitions = gft.partitions(N)
	for i in partitions[0:len(partitions)/2]:
		pylab.axvline((N/2+i)/float(N),color='r',alpha=0.2)
		pylab.axvline((N/2-i)/float(N),color='r',alpha=0.2)


# Plot a 2D image and it's GFT, with windows and partitions
def plot2D(image,window,figTitle=None):
	N = numpy.shape(image)[0]
	transform = gft.gft2d(image,window)
	pylab.figure()
	pylab.subplot(121)
	if figTitle: pylab.title(figTitle)
	pylab.imshow(image.real)

	pylab.subplot(122)
	pylab.imshow(abs(transform))
	partitions = gft.partitions(N)
	for i in partitions[0:len(partitions)/2]:
		pylab.axvline((N/2+i),color='r',alpha=0.2)
		pylab.axvline((N/2-i),color='r',alpha=0.2)
		pylab.axhline((N/2+i),color='r',alpha=0.2)
		pylab.axhline((N/2-i),color='r',alpha=0.2)
	pylab.axis((0,N,0,N))



# Begin examples


# 1D tests
def transform1D():
	N = 256
	x = numpy.arange(0,1,1./N)
	signal = numpy.zeros(len(x),numpy.complex128)

	signal1 = numpy.copy(signal); signal1[N/2] = 1.0
	plot1D(signal1,'gaussian','Transform of a Delta: Gaussian Window')
	plot1D(signal1,'box','Transform of a Delta: Boxcar Window')
	
	signal2 = numpy.copy(signal); signal2[N/4] = 1.0; signal2[N/4+N/2] = 1.0
	plot1D(signal2,'gaussian','Transform of Two Deltas: Gaussian Window')
	plot1D(signal2,'box','Transform of Two Deltas: Boxcar Window')
	
	signal3 = numpy.sin(2*numpy.pi*24*x)*gaussianWindow(signal,N/4,8) + numpy.sin(2*numpy.pi*96*x)*gaussianWindow(signal,3*N/4,8) + numpy.sin(2*numpy.pi*80*x)*gaussianWindow(signal,3*N/4,8)
	plot1D(signal3,'gaussian','Transform of Three Modulated Sines: Gaussian Window')
	plot1D(signal3,'box','Transform of Three Modulated Sines: Boxcar Window')

	signal4 = numpy.sin(2*numpy.pi*(x)*(x)*N/4)
	plot1D(signal4,'gaussian','Transform of A Sine of Increasing Frequency: Gaussian Window')
	plot1D(signal4,'box','Transform of A Sine of Increasing Frequency: Boxcar Window')
	
	
# 2D tests
def transform2D():
	N = 256
	image = numpy.zeros((N,N))

	image1 = numpy.copy(image); image1[N/2,N/2] = 1.0
	plot2D(image1,'gaussian','Transform of A 2D Delta: Gaussian Window')
	plot2D(image1,'box','Transform of A 2D Delta: Boxcar Window')

	image2 = dist(image); image2 /= numpy.max(numpy.ravel(numpy.abs(image2))); image2 = numpy.cos(2*numpy.pi*96*image2)
	plot2D(image2,'gaussian','Transform of A Circularly Symmetric Cosine: Gaussian Window')
	plot2D(image2,'box','Transform of A Circularly Symmetric Cosine: Boxcar Window')

	image3 = dist(image); image3 /= numpy.max(numpy.ravel(numpy.abs(image3))); image3 = numpy.cos(2*numpy.pi*48*image3)*gaussian2D(image3,8)
	image3 = gft.shift(image3,-N/4,-N/4)
	image3 += gft.shift(image2*gaussian2D(image3,8),N/4,N/4)
	plot2D(image3,'gaussian','Transform of A Gaussian Modulated Circularly Symmetric Cosine: Gaussian Window')
	plot2D(image3,'box','Transform of A Gaussian Modulated Circularly Symmetric Cosine: Boxcar Window')


# Run the examples
def run():
	transform1D()
	transform2D()
	
	
# simple example
def simpleExample():
	# Create a signal
	N = 256
	signal = numpy.zeros(N,numpy.complex128)
	# Make our signal 1 in the middle, zero elsewhere
	signal[N/2] = 1.0
	# plot our signal
	pylab.figure()
	pylab.plot(abs(signal))
	# Do the GFT on the signal
	SIGNAL = gft.gft1d(signal,'gaussian')
	# plot the magnitude of the signal
	pylab.figure()
	pylab.plot(abs(SIGNAL),'b')
	
	# get the partitions
	partitions = gft.partitions(N)
	# for each partition, draw a line on the graph.  Since the partitions only indicate the
	# positive frequencies, we need to draw partitions on both sides of the DC
	for i in partitions[0:len(partitions)/2]:
		pylab.axvline((N/2+i),color='r',alpha=0.2)
		pylab.axvline((N/2-i),color='r',alpha=0.2)
		
	# finally, interpolate the GFT spectrum and plot a spectrogram
	pylab.figure()
	pylab.imshow(abs(gft.gft1dInterpolateNN(SIGNAL)))
	
	
	
	
# if we're being run from the command line, execute the run function
if (__name__ == '__main__'):
	run()

