/** -----------------------------------------------------------
  MEX file for computing a solution to the LSAP with the Hungarian algorithm
 
             [rho] = hungarianLSAP(C,init_type,forb)
      [rho,varrho] = hungarianLSAP(C,init_type,forb)
         [rho,u,v] = hungarianLSAP(C,init_type,forb)
  [rho,varrho,u,v] = hungarianLSAP(C,init_type,forb)

  Given a nxm cost matrix C (integer ou floatting values) it computes:

  - a solution rho to the LSAP, i.e. a one-to-one mapping from the rows of C to its columns,
    and the one-to-to mapping from the columns to the rows

    rho is represented by a nx1 matrix so that rho(i)=j means that i is assigned to j
    varrho is a 1xm matrix so that varrho(j)=i means that j is assigned to i
    for the asymmetric LSAP (n and m are not equal):
      if n<m:
         rho is an injection from {1,...,n} to {1,...,m}
         m-n columns are not assigned, which is represented by varrho(j)=-1
      if n>m
         varrho is an injection from {1,...,m} to {1,...,n}
         n-m rows are not assigned, which is represented by rho(i)=-1

  - a solution (u,v) to its dual problem (labeling problem)

  optional integer init_type:
  0 -> no initialization
  1 -> classical initialization u is the min of each row and v the min of each column minus u

  optional boolean parameter forb:
  - true  -> forbidden assignments are represented by negative cost values
  - false -> no forbidden assignments (by default) 

  This is a MEX-file for MATLAB.
  
  This file is part of LSAPE.
  LSAPE is free software: you can redistribute it and/or modify it
  under the terms of the CeCILL-C License. See README for more details.

     Copyright 2015-2017
      authors: Sebastien Bougleux
  institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, France
   last modif: July 5 2017
 
  execute matlab file 'compile_mex' to compile this function and use it in matlab
 *  -----------------------------------------------------------
*/

#include <mex.h>
#include <cstdio>
#include <string>
#include <hungarian-lsap.hh>
#include <utils.hh>

template <typename DT>
void hungarian_helper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const mxClassID &mxC)
{
  // 1st input : edit cost matrix
  const DT *C = (DT*)mxGetPr(prhs[0]);
  // get dimensions
  int nrows = (int)mxGetDimensions(prhs[0])[0];
  int ncols = (int)mxGetDimensions(prhs[0])[1];
  int dimr2[2] = { nrows, 1 }, dimc2[2] = { 1, ncols };
  bool forbassign = false;
  unsigned short init_type = 1;
  if (nrhs > 1) init_type = (unsigned short)mxGetScalar(prhs[1]);
  if (nrhs == 3) forbassign = (bool)mxGetScalar(prhs[2]);
  //------------------------------------------------------------------
  int i = 0;
  plhs[0] = mxCreateNumericArray(2, dimr2, mxINT32_CLASS, mxREAL);
  int *rho = (int*)mxGetData(plhs[0]);
  int *varrho = NULL;
  if (nlhs == 2 || nlhs == 4)
  {
    plhs[1] = mxCreateNumericArray(2, dimc2, mxINT32_CLASS, mxREAL);
    varrho = (int*)mxGetData(plhs[1]);
  }
  DT *u = NULL, *v = NULL;
  if (nlhs > 2)
  {
    plhs[nlhs-2] = mxCreateNumericArray(2, dimr2, mxC, mxREAL);
    u = (DT*)mxGetPr(plhs[nlhs-2]);
    plhs[nlhs-1] = mxCreateNumericArray(2, dimc2, mxC, mxREAL);
    v = (DT*)mxGetPr(plhs[nlhs-1]);
  }
  else { u = new DT[nrows]; v = new DT[ncols];}
  hungarianLSAP<DT,int>(C,nrows,ncols,rho,u,v,varrho,init_type,forbassign);
  for (int i = 0; i < nrows; i++) rho[i]++;
  if (nlhs == 2 || nlhs == 4) for (int j = 0; j < ncols; j++) varrho[j]++;
  if (nlhs < 3) { delete[] u; delete[] v; }
}

//------------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  if (nrhs < 1 || nrhs > 3 || nlhs < 1 || nlhs > 4) mexErrMsgTxt("USAGE: [rho] = hungarianLSAP(C,init_type,forbassign)\nUSAGE: [rho,varrho] = hungarianLSAP(C,init_type,forbassign)\nUSAGE: [rho,u,v] = hungarianLSAP(C,init_type,forbassign)\nUSAGE: [rho,varrho,u,v] = hungarianLSAP(C,init_type,forbassign)\n");
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
