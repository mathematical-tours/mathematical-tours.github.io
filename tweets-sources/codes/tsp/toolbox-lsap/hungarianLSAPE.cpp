/** -----------------------------------------------------------
 * Matlab interface for solving the LSAPE
 *
 * -----------------------------------------------------------
 * authors: Sebastien Bougleux
 * institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France
 *
 * ----------------------------------------------------------- 
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * ----------------------------------------------------------- 
 *
 *         [rho,varrho] = hungarianLSAPE(C,init_type,forbassign)
 *     [rho,varrho,u,v] = hungarianLSAPE(C,init_type,forbassign)
 * [rho,varrho,mincost] = hungarianLSAPE(C,init_type,forbassign)
 * 
 * Given an edit cost matrix C, compute a solution (rho,varrho) to the LSAPE 
 * and a solution (u,v) to its dual problem
 * init_type: 0 (no initialization), 1 (basic initialization)
 * forbassign is an optional boolean (false by default):
 *   - true -> some value are negative in C (forbidden assignments)
 *   - false -> no forbidden assignments
 *
 * -----------------------------------------------------------
 * execute matlab file 'compile_mex' to compile this function
 * and use it in matlab
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <cstdio>
#include <string>
#include <hungarian-lsape.hh>

template <typename DT>
void hungarian_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  // get dimensions
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  int nr1 = nrows-1, nc1 = ncols-1;
  int dimr2[2] = { nrows, 1 }, dimc2[2] = { 1, ncols };
  bool forbassign = false;
  unsigned short init_type = 1;
  if (nrhs >= 2) init_type = (unsigned short)mxGetScalar(prhs[1]);
  if (nrhs == 3) forbassign = (bool)mxGetScalar(prhs[2]);
  //------------------------------------------------------------------
  int dimr[2] = { nr1, 1 }, dimc[2] = { 1, nc1 }, i = 0;
  plhs[0] = mxCreateNumericArray(2, dimr, mxINT32_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(2, dimc, mxINT32_CLASS, mxREAL);
  int *rho = (int*)mxGetData(plhs[0]);
  int *varrho = (int*)mxGetData(plhs[1]);
  DT *u = NULL, *v = NULL;
  if (nlhs == 4)
  {
    plhs[2] = mxCreateNumericArray(2, dimr2, mxC, mxREAL);
    u = (DT*)mxGetPr(plhs[2]);
    plhs[3] = mxCreateNumericArray(2, dimc2, mxC, mxREAL);
    v = (DT*)mxGetPr(plhs[3]);
  }
  else { u = new DT[nrows]; v = new DT[ncols];}
  hungarianLSAPE<DT,int>(C,nrows,ncols,rho,varrho,u,v,init_type,forbassign);
  if (nlhs == 3)
  {
    DT mincost = 0;
    for (int i = 0; i < nr1; i++) { rho[i]++; mincost += u[i]; }
    for (int j = 0; j < nc1; j++) { varrho[j]++; mincost += v[j]; }
    plhs[2] = mxCreateDoubleScalar((double)mincost);
  }
  else
  {
    for (int i = 0; i < nr1; i++) rho[i]++;
    for (int j = 0; j < nc1; j++) varrho[j]++;
  }
  if (nlhs != 4) { delete[] u; delete[] v; }
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) 
{
  if (nrhs < 1 || nrhs > 3 || nlhs < 2 || nlhs > 4) mexErrMsgTxt("USAGE: [rho,varrho] = hungarianLSAPE(C,init_type,forbassign)\nUSAGE: [rho,varrho,mincost] = hungarianLSAPE(C,init_type,forbassign)\nUSAGE: [rho,varrho,u,v] = hungarianLSAPE(C,init_type,forbassign)\n");
  mxClassID tclass = mxGetClassID(prhs[0]);
  switch(tclass)
  {
    case mxINT16_CLASS: hungarian_helper<short>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT32_CLASS: hungarian_helper<int>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxINT64_CLASS: hungarian_helper<int64_T>(nlhs,plhs,nrhs,prhs,tclass); break;
    case mxSINGLE_CLASS: hungarian_helper<float>(nlhs,plhs,nrhs,prhs,tclass); break;
    default: hungarian_helper<double>(nlhs,plhs,nrhs,prhs,tclass);
  }
}
