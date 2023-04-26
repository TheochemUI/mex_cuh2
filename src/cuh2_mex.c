/*
 * cuh2_mex.c - an interface to the CuH2 potential of the Jonsson group
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

#include "cuh2_mex.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* Ensure valid inputs */
  if (nrhs != 3) {
    mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidNumInputs",
                      "Three input arguments required.");
  } else if (nlhs > 1) {
    mexErrMsgIdAndTxt("MATLAB:cuh2pot:maxlhs", "Too many output arguments.");
  }

  if (!mxIsDouble(R_IN)) {
    mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT",
                      "First input, positions, must be a real matrix.");
  } else {
    mwSize numCols = mxGetN(R_IN);
    if (!(numCols == 3)) {
      mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT",
                        "Positions must be one x,y,z per row.");
    }
  }

  if (!(mxIsInt32(ATMNRS_IN))) {
    mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT",
                      "Second input, positions, must be an integer vector.");
  } else {
    mwSize numElems = mxGetNumberOfElements(ATMNRS_IN);
    if (!(numElems == mxGetM(R_IN))) {
      mexErrMsgIdAndTxt(
          "MATLAB:cuh2pot:invalidT",
          "Atomic number vector not the same size as the number of atoms.");
    }
  }

  if (!mxIsDouble(BOX_IN)) {
    mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT",
                      "Third input, box, must be a real 3x3 matrix.");
  } else {
    mwSize numRows = mxGetM(BOX_IN);
    mwSize numCols = mxGetN(BOX_IN);
    if (numCols != 3 || numRows != 3) {
      mexErrMsgIdAndTxt("MATLAB:cuh2pot:invalidT",
                        "Third input, box, must be a real 3x3 matrix.");
    }
  }

  /* Marshall call */
  /* XXX: This is wrong, the caller needs to fix the order for the matrices */
  double *R = mxGetPr(R_IN);
  int *atomicNrs = (int *)mxGetPr(ATMNRS_IN);
  double *box_in = mxGetPr(BOX_IN);
  int natoms = mxGetM(R_IN);
  int ndim = natoms * 3;
  mxArray *forces = mxCreateDoubleMatrix(natoms, 3, mxREAL);
  double *F = mxGetPr(forces);
  double e_val = 0;
  // Count atoms
  int *natms = (int *)malloc(2 * sizeof(int)); // Always Cu, then H
  natms[0] = 0;
  natms[1] = 0;
  // Simple counter
  for (int idx = 0; idx < natoms; idx++) {
    if (atomicNrs[idx] == 29) {
      natms[0] += 1; // Cu
    } else if (atomicNrs[idx] == 1) {
      natms[1] += 1; // H
    } else {
      mexErrMsgIdAndTxt(
          "MATLAB:cuh2pot:invalidT",
          "Only 29(Cu) and 1 (H), are allowed in the atomic numbers input");
    }
  }

  /* We need to transpose the inputs we got */
  /* Transpose R */
  double *R_transposed =
      (double *)malloc(mxGetN(R_IN) * mxGetM(R_IN) * sizeof(double));
  for (int idx = 0; idx < mxGetN(R_IN); idx++) {
    for (int jdx = 0; jdx < mxGetM(R_IN); jdx++) {
      R_transposed[jdx * mxGetN(R_IN) + idx] = R[idx * mxGetM(R_IN) + jdx];
    }
  }
  /* Only diagonal boxes are supported for CuH2 */
  double box[3] = {box_in[0], box_in[4], box_in[8]};

  /* Call! */
  c_force_eam(natms, ndim, box, R_transposed, F, &e_val);

#ifdef DEBUG
  /* Print values */
  char buf[100];
  sprintf(buf, "natoms = %d, ndim = %d\nbox = %f %f %f\n", natoms, ndim, box[0],
          box[1], box[2]);
  printf("%s", buf);
  printf("F:\n");
  for (int idx = 0; idx < natoms; idx++) {
    for (jdx = 0; jdx < 3; jdx++) {
      printf("%f ", F[idx * 3 + jdx]);
    }
    printf("\n");
  }

  printf("R:\n");
  for (int idx = 0; idx < natoms; idx++) {
    for (jdx = 0; jdx < 3; jdx++) {
      printf("%f ", R_transposed[idx * 3 + jdx]);
    }
    printf("\n");
  }

  printf("e_val: %f\n", e_val);
#endif

  /* MATLAB structure for outputs */
  const char *onames[] = {"energy", "forces"};
  mxArray *oarray = mxCreateStructMatrix(1, 1, 2, onames);
  // set energy to a double
  mxArray *energy = mxCreateDoubleScalar(e_val);
  mxSetField(oarray, 0, "energy", energy);

  /* set forces to a real matrix */
  mxSetField(oarray, 0, "forces", forces);

  EF_OUT = oarray;
}
