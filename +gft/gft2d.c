#include <string.h>
#include "mex.h"
#include "gft.h"

/* gft2d(sig, windowType)
 *
 * Performs a 2d GFT
 * (doesn't shift the output)
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1 || nrhs > 2) 
        mexErrMsgTxt("gft2d: Wrong number of inputs");
    else if (nlhs>1)
        mexErrMsgTxt("gft2d: Too many output arguments");

    int winType = 0; /* 0=gaussian, 1=box */
    if (nrhs==2) 
    {
        /* user specified window type */
        if (!mxIsChar(prhs[1]))
            mexErrMsgTxt("gft2d: windowType must be a string");
        
        char *windowType;
        windowType = mxArrayToString(prhs[1]);
        if (windowType==NULL)
            mexErrMsgTxt("gft2d: problem with windowType");
        if (strcmp(windowType, "gaussian"))
            winType = 1;
    }

    mwSize n = mxGetNumberOfDimensions(prhs[0]);
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || n!=2)
        mexErrMsgTxt("gft2d: only real double 2d arrays supported");
    mwSize m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    mwSize Nele = m*n;

    double *inpr = mxGetPr(prhs[0]);

    int *pars = gft_1dPartitions(Nele);
    void *window;
    if (winType==0)
        window = &gaussian;
    else if (winType==1)
        window = &box;
    else
        mexErrMsgTxt("gft2d: unknown window type");
        
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
    /* tranpose n/m here since Matlab has Fortran order */
    gft_2dComplex64(cdata, n, m, window);
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

    free(cdata);
    plhs[0] = out;

}
