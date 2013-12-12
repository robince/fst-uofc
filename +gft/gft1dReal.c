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
        mexErrMsgTxt("gft1dreal: Wrong number of inputs");
    else if (nlhs>1)
        mexErrMsgTxt("gft1dreal: Too many output arguments");

    mwSize n = mxGetNumberOfDimensions(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || n!=2)
        mexErrMsgTxt("gft1dreal: only real double 2d arrays supported");
    mwSize m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (MIN(m,n) != 1)
        mexErrMsgTxt("gft1d: input should be a vector");

    mwSize Nele = mxGetNumberOfElements(prhs[0]);
    double *inpr = mxGetPr(prhs[0]);

    int *pars = gft_1dRealPartitions(Nele);
    double *win = windowsFromPars(Nele, &gaussian, pars);
        
    /* perform computation
     * cdata is complex interlaced format 
     */
    double *cdata;
    cdata = (double *)malloc(Nele*sizeof(double)*2);
    mwSize i;
    for(i=0; i<Nele; i++)
    {
        cdata[(2*i)] = inpr[i]; /* real part */
        cdata[(2*i)+1] = 0.0;   /* imag part */
    }
    /* do GFT computation in place */
    gft_1dComplex64(cdata, Nele, win, pars, 1);
    /* shift */
    /* shift(cdata, Nele, Nele/2); */

    mxArray *out = mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxCOMPLEX);
    double *outpr = mxGetPr(out);
    double *outpi = mxGetPi(out);

    for(i=0; i<Nele; i++)
    {
        outpr[i] = cdata[(2*i)];   /* real */
        outpi[i] = cdata[(2*i)+1]; /* imag */
    }

    free(win);
    free(pars);
    free(cdata);
    plhs[0] = out;

}
