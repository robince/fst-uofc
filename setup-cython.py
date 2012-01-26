#	This software is copyright (c) 2010 UTI Limited Partnership.  
#	The original authors are Robert A. Brown, M. Louis Lauzon 
#	and Richard Frayne.  This software is licensed in the terms 
#	set forth in the "FST License Notice.txt" file, which is 
#	included in the LICENSE directory of this distribution.


from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


sourcefiles = ['pyGFT.pyx','gft.c']

ext_modules = [
	Extension("gft", sourcefiles,
						libraries=['fftw3'],
						include_dirs = ['/usr/local/include']
						)
	]

setup(
  name = 'gft',
	version = "0.1",
	description="A Python module to calculate the General Fourier Family Transform",
	author_email="robb@robbtech.com",
  cmdclass = {'build_ext': build_ext},
	include_dirs = [numpy.get_include()],
  ext_modules = ext_modules
)
