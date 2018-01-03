#include <iostream>
#include <math.h>
#include <map>
#include <omp.h>
#include <cstring>
#include <vector>
#include <cassert>
#include "mex.h"
using namespace std;

// Adapted from Andrew Adams' ImageStack library

extern void _main();

void calculateCoefficients(double sigma, double *c0, double *c1, double *c2, double *c3) {
    // performs the necessary conversion between the sigma of a Gaussian blur
    // and the coefficients used in the IIR filter

    double q;

    assert(sigma >= 0.5);

    if (sigma < 2.5) {
        q = (3.97156 - 4.14554*sqrt(1 - 0.26891*sigma));
    } else {
        q = 0.98711*sigma - 0.96330;
    }

    double denom = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    *c1 = (2.44413*q + 2.85619*q*q + 1.26661*q*q*q)/denom;
    *c2 = -(1.4281*q*q + 1.26661*q*q*q)/denom;
    *c3 = (0.422205*q*q*q)/denom;
    *c0 = 1 - (*c1 + *c2 + *c3);
}

void computeAttenuation(double *scale, int size, int width,
                                  const double c0, const double c1,
                                  const double c2, const double c3,
                                  int iterations) {
    // Initial value
    for (int x = 0; x < width; x++) scale[x] = 1.0;
    for (int x = width; x < size; x++) scale[x] = 0.0;

    for (int i = 0; i < iterations; i++) {
        // Forward pass
        scale[0] = c0 * scale[0];
        scale[1] = c0 * scale[1] + c1 * scale[0];
        scale[2] = c0 * scale[2] + c1 * scale[1] + c2 * scale[0];
        for (int x = 3; x < size; x++) {
            scale[x] = (c0 * scale[x] +
                        c1 * scale[x-1] +
                        c2 * scale[x-2] +
                        c3 * scale[x-3]);
        }

        // Backward pass
        scale[size-1] = c0 * scale[size-1];
        scale[size-2] = (c0 * scale[size-2] +
                         c1 * scale[size-1]);
        scale[size-3] = (c0 * scale[size-3] +
                         c1 * scale[size-2] +
                         c2 * scale[size-1]);
        for (int x = size-4; x >= 0; x--) {
            scale[x] = (c0 * scale[x] +
                        c1 * scale[x+1] +
                        c2 * scale[x+2] +
                        c3 * scale[x+3]);
        }
    }

    // Invert
    for (int x = 0; x < size; x++) scale[x] = 1.0/scale[x];
}

void blurChunk(double *data, int size,
                         const double c0, const double c1,
                         const double c2, const double c3) {

    data[0] *= c0;
    data[1] = c0*data[1]+c1*data[0];
    data[2] = c0*data[2]+c1*data[1]+c2*data[0];

    // Filter
    for (int i = 3; i < size; i++) {
        data[i] = (c0 * data[i] +
                   c1 * data[i-1] +
                   c2 * data[i-2] +
                   c3 * data[i-3]);
    }
    
    data[size-1] *= c0;
    data[size-2] = c0*data[size-2]+c1*data[size-1];
    data[size-3] = c0*data[size-3]+c1*data[size-2]+c2*data[size-1];
    for (int i = size-4; i >= 0; i--) {
        data[i] = (c0 * data[i] +
                   c1 * data[i+1] +
                   c2 * data[i+2] +
                   c3 * data[i+3]);
    }
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray *prhs[]) {
	if (nrhs != 2)
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", "Requires image and standard deviation.");
    else if (nlhs != 1) 
        mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout", "Produces an image.");

    // Read list of vertices
    double *image = mxGetPr(prhs[0]);
    int r = (int)mxGetM(prhs[0]);
    int c = (int)mxGetN(prhs[0]);
    
    double filterWidth = *(double*)mxGetPr(prhs[1]);
    double *result = (double*)mxCalloc(r*c,sizeof(double));
    
    // DO MATH ////////////////////////////////////////////////////////////
    
    int iterations = 1;
    while (filterWidth > 64) {
        filterWidth /= sqrtf(2);
        iterations *= 2;
    }
    
    const int size = r + (int)(filterWidth*6);

    double c0, c1, c2, c3;
    calculateCoefficients(filterWidth, &c0, &c1, &c2, &c3);

    vector<double> scale(size);
    computeAttenuation(&scale[0], size, r, c0, c1, c2, c3, iterations);

    #pragma omp parallel for
    for (int col = 0; col < c; col++) {
        vector<double> chunk(size);
        memcpy(&chunk[0], &image[col*r], r*sizeof(double));
        
        for (int i = 0; i < iterations; i++)
            blurChunk(&chunk[0], size, c0, c1, c2, c3);
        
        for (int row = 0; row < r; row++)
            result[row+col*r] = chunk[row]*scale[row];
    }

    ///////////////////////////////////////////////////////////////////////

    plhs[0] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
    mxSetData(plhs[0],result);
    mxSetM(plhs[0],r);
    mxSetN(plhs[0],c);
}