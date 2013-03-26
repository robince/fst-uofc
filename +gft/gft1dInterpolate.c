#include <string.h>
#include <math.h>
#include "mex.h"
#include "gft.h"

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/* gft1dInterpolate(st, M)
 *
 * Interpolate a 1d GFT to show spectrogram
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1 || nrhs > 2) 
        mexErrMsgTxt("gft1d: Wrong number of inputs");
    else if (nlhs>1)
        mexErrMsgTxt("gft1d: Too many output arguments");

    mwSize n = mxGetNumberOfDimensions(prhs[0]);
    if (!mxIsDouble(prhs[0]) || !mxIsComplex(prhs[0]) || n!=2)
        mexErrMsgTxt("gft1d: only complex double 2d arrays supported");
    mwSize m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (MIN(m,n) != 1)
        mexErrMsgTxt("gft1d: input should be a vector");
    mwSize Nele = mxGetNumberOfElements(prhs[0]);
    mexPrintf("Nele: %d \n", Nele);

    if (nrhs==2) 
    {
        /* user specified M */
        mwSize n = mxGetNumberOfElements(prhs[1]);
        if (!mxIsDouble(prhs[1]) || n != 1)
            mexErrMsgTxt("gft2d: M must be a scalar value");
        double *pr = mxGetPr(prhs[1]);
        m = round(pr[0]);
    }
    else
    {
        m = Nele;
    }
    
    double *inpr = mxGetPr(prhs[0]);
    double *inpi = mxGetPi(prhs[0]);
        
     /* cdata is complex interlaced format 
     */
    double *cdata;
    cdata = (double *)malloc(Nele*sizeof(double)*2);
    mwSize i;
    for(i=0; i<Nele; i++)
    {
        cdata[(2*i)] = inpr[i]; /* real part */
        cdata[(2*i)+1] = inpi[i];   /* imag part */
    }

    /* shift */
    shift(cdata, Nele, -Nele/2);

    /* perform interpolation */
    unsigned int Nint = Nele;
    unsigned int Mint = m;
    double *image = gft_1d_interpolateNN(cdata, Nint, 0);

    mxArray *out = mxCreateNumericMatrix(m, m, mxDOUBLE_CLASS, mxCOMPLEX);
    double *outpr = mxGetPr(out);
    double *outpi = mxGetPi(out);

    for(i=0; i<(m*m); i++)
    {
        outpr[i] = image[(2*i)];   /* real */
        outpi[i] = image[(2*i)+1]; /* imag */
        /* debugging problem with NaNs and high values in output */
        if(image[(2*i)]>1000)
            mexPrintf("loop: %d %f ... %f\n",i,image[(2*i)],outpr[i]);
    }

    free(cdata);
    free(image);
    plhs[0] = out;

}
