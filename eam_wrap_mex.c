/*
 * eam_wrap_mex.c - an interface to the CuH2 potential of the Jonsson group
 *
 *  Created on: 25 April 2023
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 *     License: MIT
 *
 * Multiplies an input scalar (multiplier)
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *    [E, F] = cuh2pot(posMatrix, atmNrsMat, boxMat)
 *
 * This is a MEX file for MATLAB.
 */

#include "eam_wrap_mex.h"

#define EF_OUT plhs[0]

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *R, *F, *U;
    size_t *atomicNrs;
    size_t m, n, natoms;

    // Ensure valid inputs
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidNumInputs", "Three input arguments required.");
    } else if (nlhs > 2) {
        mexErrMsgIdAndTxt("MATLAB:cuh2pot:maxlhs", "Too many output arguments.");
    }

    /* if (!mxIsDouble(R_IN)) { */
    /*     mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT", "First input, positions, must be a real matrix."); */
    /* } */

    /* if (!( mxIsInt32(ATMNRS_IN) )) { */
    /*     mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT", "Second input, positions, must be an integer vector."); */
    /* } */

    if (!mxIsDouble(BOX_IN)) {
        mwSize numRows = mxGetM(BOX_IN);
        mwSize numCols = mxGetN(BOX_IN);
        if (!(numRows == numCols == 3)){
            mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT", "Third input, box, must be a real 3x3 matrix.");
        }
    }

    // MATLAB structure for outputs
    const char *onames[] = {"energy", "forces"};
    mxArray *oarray = mxCreateStructMatrix(1, 1, 2, onames);
    // set energy to a double
    double val = 3.14;
    mxArray *energy = mxCreateDoubleScalar(val);
    mxSetField(oarray, 0, "energy", energy);

    // set forces to a real matrix
    mxArray *forces = mxCreateDoubleMatrix(1, 3, mxREAL);
    double *data = mxGetPr(forces);
    data[0] = 1.0;
    data[1] = 3.0;
    data[2] = 7.0;
    mxSetField(oarray, 0, "forces", forces);

    EF_OUT = oarray;
}
