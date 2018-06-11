// =========================================================================
/** \file greedy-lsap.hh
 *  \brief Greedy algorithms for approximating symmetric and asymmetric Linear Sum Assignment Problems (LSAP)
 * \author Sebastien Bougleux (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France)
 */
// =========================================================================
/* 
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.
   
   You should have received a copy of the CeCILL-C License along with 
   LSAPE. If not, see http://www.cecill.info/index.en.html

   -----------------------------------------------------------
   
   Creation: December 5 2015
   Modifications: July 6 2017
   
   -----------------------------------------------------------
   Main function:

   cost = greedyLSAP(C,n,m,rho,varrho=NULL,greedy_type=2)

   -> see at the end of this file or generate doxygen documentation
*/
// =========================================================================

#ifndef _GREEDY_LSAP_
#define _GREEDY_LSAP_

#include <iostream>
#include <limits>
#include <algorithm>
#include "utils.hh"

// -----------------------------------------------------------
// Basic greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyBasicLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{
  IT i = -1, imin, itmp;
  DT cmin, mxdt = std::numeric_limits<DT>::max(), approx = 0;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }

  IT *unass = new IT[n+1], *pti_unass = NULL;
  IT *pti_unass_beg = unass, *pti_unass_end = unass+n, *pti_min = NULL;
  
  for (i = 0; i < n; i++) { unass[i] = i; rho[i] = -1; }
  
  // augmentation of columns
  for (IT j = 0; j < m; j++)
  {
    // find the min among unassigned rows
    cmin = mxdt;
    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
    {
      const DT &cij = C[j*n+*pti_unass];
      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
    }
    // assign the row which provides the minimum and update the approximate solution
    imin = *pti_min; rho[imin] = j; varrho[j] = imin;
    *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
    approx += cmin;
  }
  
  delete[] unass;
  if (deletevarrho) delete[] varrho;
  return approx;
}

// -----------------------------------------------------------
// Refined greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyRefinedLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{
  IT nass = 0, i = -1, j = -1, imin, jmin;
  DT cmin, ckmin, mxdt = std::numeric_limits<DT>::max(), approx = 0;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }

  IT *unassi = new IT[n+1], *pti_unass = NULL;
  IT *pti_unass_beg = unassi, *pti_unass_end = unassi+n, *pti_min = NULL;

  IT *ptj_unass = NULL, *unassj = new IT[m+1];
  IT *ptj_unass_beg = unassj, *ptj_unass_end = unassj+m, *ptj_min = NULL;

  IT *ptj_unass1 = NULL;
  
  for (i = 0; i < n; i++) { unassi[i] = i; rho[i] = -1; }
  for (j = 0; j < m; j++) { unassj[j] = j; }
  
  // augmentation of columns
  for (ptj_unass1 = ptj_unass_beg; ptj_unass1 != ptj_unass_end;)
  {
    j = *ptj_unass1;
    // find the min among unassigned rows
    cmin = mxdt;
    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
    {
      const DT &cij = C[j*n+*pti_unass];
      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
    }
    // find the min among unassigned columns for imin
    imin = *pti_min;
    ckmin = mxdt;
    for (ptj_unass = ptj_unass_beg; ptj_unass != ptj_unass_end; ptj_unass++)
    {
      const DT &cik = C[*ptj_unass*n+imin];
      if (cik  < ckmin) { ckmin = cik; ptj_min = ptj_unass; }
    }
    // assign the row and column which provides the minimum
    if (cmin <= ckmin)
    {
      rho[imin] = j; varrho[j] = imin;
      *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
      ptj_unass_beg++; 
      approx += cmin;
    }
    else
    {
      jmin = *ptj_min; rho[imin] = jmin; varrho[jmin] = imin;
      *ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
      *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
      approx += ckmin;
    }
    ptj_unass1 = ptj_unass_beg;
  }
  
  delete[] unassi; delete[] unassj;
  if (deletevarrho) delete[] varrho;
  return approx;
}

// -----------------------------------------------------------
// Loss greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyLossLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{
  IT nass = 0, i = -1, j = -1, imin, jmin, imin2, imin3;
  DT cmin, cij, ckmin, mxdt = std::numeric_limits<DT>::max(), approx = 0, cmin2, cmin3;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }

  IT *unassi = new IT[n+1], *pti_unass = NULL;
  IT *pti_unass_beg = unassi, *pti_unass_end = unassi+n, *pti_min = NULL, *pti_min2 = NULL, *pti_min3 = NULL;

  IT *ptj_unass = NULL, *unassj = new IT[m+1];
  IT *ptj_unass_beg = unassj, *ptj_unass_end = unassj+m, *ptj_min = NULL;

  IT *ptj_unass1 = NULL;
  
  for (i = 0; i < n; i++) { unassi[i] = i; rho[i] = -1; }
  for (j = 0; j < m; j++) { unassj[j] = j; }
  
  // augmentation of columns
  for (ptj_unass1 = ptj_unass_beg; ptj_unass1 != ptj_unass_end;)
  {
    j = *ptj_unass1;
    // find the min among unassigned rows
    cmin = mxdt;
    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
    {
      cij = C[j*n+*pti_unass];
      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
    }
    // find the min among unassigned columns for imin
    imin = *pti_min;
    ckmin = mxdt;
    for (ptj_unass = ptj_unass_beg; ptj_unass != ptj_unass_end; ptj_unass++)
    {
      const DT &cik = C[*ptj_unass*n+imin];
      if (cik  < ckmin) { ckmin = cik; ptj_min = ptj_unass; }
    }
    // assign the row and column which provides the minimum
    jmin = *ptj_min;
    if (j == jmin)
    {
      rho[imin] = j; varrho[j] = imin;
      *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
      ptj_unass_beg++;
      approx += cmin;
    }
    else
    {
      // find the min among unassigned rows different to imin, for j => find the 2nd min
      cmin3 = cmin2 = mxdt;
      for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
      {
	if (*pti_unass != imin)
	{
	  cij = C[j*n+*pti_unass];
	  if (cij  < cmin2) { cmin2 = cij; pti_min2 = pti_unass; }
	  cij = C[jmin*n+*pti_unass];
	  if (cij  < cmin3) { cmin3 = cij; pti_min3 = pti_unass; }
	}
      }
      imin2 = *pti_min2;
      imin3 = *pti_min3;
      if (cmin + cmin3 < cmin2 + ckmin) // remove j, jmin, imin, imin3
      {
	rho[imin] = j; varrho[j] = imin;
	rho[imin3] = jmin; varrho[jmin] = imin3;
	*pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	ptj_unass_beg++;
	*ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
	*pti_min3 = *pti_unass_beg; *pti_unass_beg = imin3; pti_unass_beg++;
	approx += cmin + cmin3;
      }
      else // remove j, jmin, imin, imin2
      {
	rho[imin2] = j; varrho[j] = imin2;
	rho[imin] = jmin; varrho[jmin] = imin;
	ptj_unass_beg++;
	*ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
	*pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	*pti_min2 = *pti_unass_beg; *pti_unass_beg = imin2; pti_unass_beg++;
	approx += ckmin + cmin2;
      }
    }
    ptj_unass1 = ptj_unass_beg;
  }
  
  delete[] unassi; delete[] unassj;
  if (deletevarrho) delete[] varrho;
  return approx;
}

// -----------------------------------------------------------
// BasicSort greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
class CostSortComp
{
private:
  const DT *_C;
public:
  CostSortComp(const DT *C) : _C(C) { }
  ~CostSortComp() { }
  inline bool operator()(const IT &ij, const IT &kl) { return (_C[ij] < _C[kl]); }
};

template <class DT, typename IT>
DT greedyBasicSortLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{ 
  // sort the costs
  IT *idxs = new IT[n*m+1];
  IT *pt_idxs_end = idxs+n*m, i = 0, j;
  DT approx = 0;
  for (IT *pt_idx = idxs; pt_idx != pt_idxs_end; pt_idx++, i++) *pt_idx = i;
  
  CostSortComp<DT,IT> comp(C);
  std::sort(idxs,idxs+n*m,comp);

  // assign element in ascending order of the costs
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }
  for (IT i = 0; i < n; i++) rho[i] = -1;
  for (IT j = 0; j < m; j++) varrho[j] = -1;
  for (IT *pt_idx = idxs, nbe = 0; pt_idx != pt_idxs_end && nbe < m; pt_idx++)
  {
    ind2sub(*pt_idx,n,i,j);
    if (rho[i] == -1 && varrho[j] == -1)
    {
      rho[i] = j;
      varrho[j] = i;
      approx += C[*pt_idx];
      nbe++;
    }
  }
  if (deletevarrho) delete[] varrho;
  delete[] idxs;
  return approx;
}

// -----------------------------------------------------------
// Counting sort greedy LSAP for n >= m only (not checked)
// and for integer values
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyCountingSortLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{
  // find min and max values
  DT minc = C[0], maxc = C[0];
  const DT *ite = C+n*m;
  for (const DT *itc = C+1; itc < ite; itc++)
  {
    const DT &vc = *itc;
    if (vc < minc) minc = vc;
    else if (vc > maxc) maxc = vc;
  }

  // construct histogram
  IT nbins = maxc - minc + 1;
  IT *bins = new IT[nbins];
  const IT *itbe = bins+nbins;
  for (IT *itb = bins; itb < itbe; itb++) *itb = 0;
  for (const DT *itc = C; itc < ite; itc++) bins[*itc-minc]++;

  // starting index for each cost value
  IT tot = 0, oldcount;
  for (IT i = 0; i < nbins; i++) { oldcount = bins[i]; bins[i] = tot; tot += oldcount; }

  // reoder the costs, preserving order of C with equal keys
  IT *idxs = new IT[n*m+1], k = 0;
  for (const DT *itc = C; itc < ite; itc++, k++) { idxs[bins[*itc-minc]] = k; bins[*itc-minc]++; }

  // assign element in ascending order of the costs
  IT *pt_idxs_end = idxs+n*m, i = 0, j, approx = 0;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }
  for (IT i = 0; i < n; i++) rho[i] = -1;
  for (IT j = 0; j < m; j++) varrho[j] = -1;
  for (IT *pt_idx = idxs, nbe = 0; pt_idx != pt_idxs_end && nbe < m; pt_idx++)
  {
    ind2sub(*pt_idx,n,i,j);
    if (rho[i] == -1 && varrho[j] == -1)
    {
      rho[i] = j;
      varrho[j] = i;
      approx += C[*pt_idx];
      nbe++;
    }
  }
  if (deletevarrho) delete[] varrho;
  delete[] idxs; delete[] bins;
  return approx;
}

template<> float greedyCountingSortLSAP<>(const float *C, const int &n, const int &m, int *rho, int *varrho) { }
template<> double greedyCountingSortLSAP<>(const double *C, const int &n, const int &m, int *rho, int *varrho) { }

// --------------------------------------------------------------------------------
// Main function: greedy algorithms for both square and rectangular cost matrices
/**
 * \brief Compute an approximate solution to the LSAP as an assignment with low cost, with greedy algorithms given a square or a rectangular cost matrix
 * \param[in] C nxm cost matrix, with n>=m, represented as an array of size \c nm obtained by concatenating its columns
 * \param[in] n Number of rows of \p C (size of the 1st set)
 * \param[in] m Number of columns of \p C (size of the 2nd set)
 * \param[out] rho Array of size \p n (must be previously allocated) for the assignment of the rows to the columns, rho[i]=-1 indicates that i is not assigned, else rho[i]=j indicates that i is assigned to j
 * \param[out] varrho Array of size \p m (must be previously allocated) for the assignement of the columns to the rows
 * \param[in] greedy_type 0:Basic, 1:Refined, 2: Loss (default), 3: Basic sort, 4: Counting sort (integers only)
 * \return Cost of the assignment or -1 if n<m
 *
 * \note Based on the reference\n
 * Approximate Graph Edit Distance in Quadratic Time. K. Riesen, M. Ferrer, H. Bunke. ACM Transactions on Computational Biology and Bioinformatics, 2015.
 */
// --------------------------------------------------------------------------------
template <class DT, typename IT>
DT greedyLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL, unsigned short greedy_type = 2)
{
  if (n < m) return -1;
  switch(greedy_type)
  {
  case 0: return greedyBasicLSAP<DT,IT>(C,n,m,rho,varrho);
  case 1: return greedyRefinedLSAP<DT,IT>(C,n,m,rho,varrho);
  case 3: return greedyBasicSortLSAP<DT,IT>(C,n,m,rho,varrho);
  case 4: return greedyCountingSortLSAP<DT,IT>(C,n,m,rho,varrho);
  default: return greedyLossLSAP<DT,IT>(C,n,m,rho,varrho);
  }
}

// -----------------------------------------------------------
#endif
