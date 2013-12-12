#include <string.h>
#include "mex.h"
#include "gft.h"

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/* gft1d(sig, windowType)
 *
 * Performs a 1d GFT using an asymmetric partitioning
 * scheme for real signals
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) 
        mexErrMsgTxt("gft1dRealPartitions: Wrong number of inputs");
    else if (nlhs>1)
        mexErrMsgTxt("gft1dRealPartitions: Too many output arguments");

    mwSize n = mxGetNumberOfDimensions(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || !mxIsScalar(prhs[0]))
        mexErrMsgTxt("gft1dRealPartitions: requires scalar argument");

    mwSize Nele;
    double *inpr = mxGetPr(prhs[0]);
	Nele = inpr[0];

    int *pars = gft_1dRealPartitions(Nele);
	int i = 0;
	int parlength;

	/* count number of elements in partition vector */
	while (pars[i]>0) i++;
	parlength = i;
			
    mxArray *out = mxCreateNumericMatrix(1, parlength, mxDOUBLE_CLASS, mxREAL);
    double *outpr = mxGetPr(out);

    for(i=0; i<parlength; i++)
    {
        outpr[i] = pars[i];
    }

    free(pars);
    plhs[0] = out;
}
