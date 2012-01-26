#	This software is copyright (c) 2010 UTI Limited Partnership.  
#	The original authors are Robert A. Brown, M. Louis Lauzon 
#	and Richard Frayne.  This software is licensed in the terms 
#	set forth in the "FST License Notice.txt" file, which is 
#	included in the LICENSE directory of this distribution.


cdef extern from "gft.h":
	int* gft_1dPartitions(unsigned int N)