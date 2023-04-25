/*
 *  Created on: 25 April 2023
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 *     License: MIT
 */
#ifndef EAM_WRAP_MEX_H_
#define EAM_WRAP_MEX_H_

#include <stdlib.h>
#include <stdio.h>

#include "mex.h"

/* Fortran ISO_C_BINDING function */
extern void c_force_eam(int *natms, int ndim, double *box, double *R, double *F,
                        double *U);

/* Entry point for MATLAB */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* MATLAB variable helpers */
typedef struct {
  double energy;
  double *forces;
} Result;

/* Input Arguments */
#define R_IN prhs[0]
#define ATMNRS_IN prhs[1]
#define BOX_IN prhs[2]
/* Output Arguments */
#define EF_OUT plhs[0]

#endif // EAM_WRAP_MEX_H_
