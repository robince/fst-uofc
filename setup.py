#	This software is copyright (c) 2010 UTI Limited Partnership.
#	The original authors are Robert A. Brown, M. Louis Lauzon
#	and Richard Frayne.  This software is licensed in the terms
#	set forth in the "FST License Notice.txt" file, which is
#	included in the LICENSE directory of this distribution.


from distutils.core import setup
from distutils.extension import Extension
import numpy


sourcefiles = ['gft.c','pyGFT.c']

ext_modules = [
	Extension("gft", sourcefiles,
						libraries=['fftw3'],
						include_dirs = ['/usr/local/include']
						)
	]

setup(
  name = 'PyGFT',
	version = "0.1",
	description="A Python module to calculate the General Fourier Family Transform",
	author_email="robb@robbtech.com",
	url='http://www.robbtech.com',
	include_dirs = [numpy.get_include()],
	ext_package = 'PyGFT',
  ext_modules = ext_modules,
	py_modules = ['PyGFT.examples']
)
