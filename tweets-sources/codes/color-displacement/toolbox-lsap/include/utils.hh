// -----------------------------------------------------------   
/** \file utils.hh
 *  \brief Manipulation of assignments/matchings
 * \author Sebastien Bougleux (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France)
*/
   
/* -----------------------------------------------------------
   
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.

   -----------------------------------------------------------
   
   Creation: December 5 2015
   Last modif: July 7 2017
   
*/

#ifndef _LSAPE_UTILS_
#define _LSAPE_UTILS_

#include <cstring>
#include <limits>

// -----------------------------------------------------------
/** \brief Transform a bipartite matching rho, from a set of size \p n to a set of size \p n, into its matrix representation X (represented as an array obtained by concatenating its columns), unassigned elements are represented in rho by negative values
 *  \param[in] rho Array of size \p n representing a partial permutation
 *  \param[in] n Number of elements in the array
 *  \param[out] X Array of size n² representing the partial permutation
 */
// -----------------------------------------------------------
template <typename IT, typename BT>
void perm2Mtx(IT *rho, const IT &n, BT *X)
{
  IT i = 0, j = 0, stop;
  BT zero = 0, one = 1;
  for (; i < n; i++)
  {
    stop = rho[i];
    if (stop < 0) for (j = 0; j < n; j++) X[j*n+i] = zero; // j not assigned
    else // j assigned
    {
      for (j = 0; j < stop; j++) X[j*n+i] = zero;
      X[j*n+i] = one;
      for (j++; j < n; j++) X[j*n+i] = zero;
    }
  }
}

// -----------------------------------------------------------
/** \brief Transform a bipartite matching rho from, a set of size \p n to a set of size \p m, into its matrix representation X (represented as an array obtained by concatenating its columns), unassigned elements are represented in rho by negative values
 *  \param[in] rho Array of size \p n representing a partial permutation
 *  \param[in] n Number of elements in the array (elements in the 1st set)
 *  \param[in] m Number of elements in the array (elements in the 2nd set)
 *  \param[out] X Array of size \p nm representing the partial permutation
*/
// -----------------------------------------------------------
template <typename IT, typename IT2, typename BT>
void pperm2Mtx(IT *rho, const IT2 &n, const IT2 &m, BT *X)
{
  IT i = 0, j = 0, stop;
  BT zero = 0, one = 1;
  for (; i < n; i++)
  {
    stop = rho[i];
    if (stop < 0) for (j = 0; j < m; j++) X[j*n+i] = zero; // j not assigned
    else // j assigned
    {
      for (j = 0; j < stop; j++) X[j*n+i] = zero;
      X[j*n+i] = one;
      for (j++; j < m; j++) X[j*n+i] = zero;
    }
  }
}

// -----------------------------------------------------------
/** \brief Reconstruct the (n+1)x(m+1) binary matrix X associated to a error-correcting bipartite matching
 * \param[in] rho Assignment from the 1st set {0,...,n} to the 2nd one {0,...,m+1}
 * \param[in] n Size of \p rho
 * \param[in] varrho Assignment from the 2nd set {0,...,m} to the 1st set {0,...,n+1}
 * \param[in] m Size of \p varrho
 * \param[out] X (n+1)x(m+1) binary matrix represented as a (n+1)(m+1) array (must be already allocated)
*/
// -----------------------------------------------------------
template <typename IT, typename BT>
void ecm2Mtx(const IT *rho, const IT &n, const IT *varrho, const IT &m, BT *X)
{
  IT i, j = 0, n1 = n+1, m1 = m+1, stop;
  BT zero = 0, one = 1, *xpt = X;
  // submatrix nxm-1
  for (; j < m; j++)
  {
    stop = varrho[j];
    if (stop < 0) for (i = 0; i < n1; i++, ++xpt) *xpt = zero; // j not assigned
    else // j assigned
    {
      for (i = 0; i < stop; i++, ++xpt) *xpt = zero;
      *xpt = one;
      ++xpt;
      for (i = stop+1; i < n1; i++, ++xpt) *xpt = zero;
    }
  }
  // last column, j=m from the previous loop, and xpt corresponds to the first element
  for (i = 0; i < n; i++, ++xpt) *xpt = (rho[i] == m ? one : zero);
  *xpt = 0; // last element
}

// -----------------------------------------------------------
/** \brief Matrix reduction from a given nxm matrix C and two vectors u and v
 * \param[in] C Matrix to reduce
 * \param[in] n Number of rows
 * \param[in] m Number of columns
 * \param[in] u Array associating a value to each row of C
 * \param[in] v Array associating a value to each column of C
 * \details
 * Update \p C by C[j*n+i] - (u[i] + v[j]) for each i=1,...,n and j=1,...,m
*/
// -----------------------------------------------------------
template <class DT, typename IT>
void reduceMtx(DT *C, const IT &n, const IT &m, DT *u, DT *v)
{
  for (IT i = 0, j = 0; i < n; i++)
  {
    const DT &ui = u[i];
    for (j = 0; j < m; j++)
    {
      if (C[j*n+i] < 0) continue;  // forbidden assignments
      C[j*n+i] -= ui + v[j];
    }
  }
}

// -----------------------------------------------------------

template <typename IT>
int sub2ind(const IT &i, const IT &j, const IT &n) {return i + j*n;}

template <typename IT>
void ind2sub(const IT &ij, const IT &n, IT &i, IT &j) { j = ij / n; i = ij - j * n; }

// -----------------------------------------------------------
/** \brief Extend a (n+1)x(m+1) LSAPE instance to an equivalent (n+m)x(m+n) LSAP instance
 *  \param[in] C nxm edit cost matrix
 *  \param[in] n number of rows (last row correspond to the null element)
 *  \param[in] m number of colums (last column correspond to the null element)
 *  \param[out] Cext (n-1+m-1)x(m-1+n-1) extended cost matrix (must be previously allocated)
 *  \param[in] valmx value given to outer diagonal elements of Cext for removals and insertions
 */
template <typename DT,typename IT>
void extendLSAPEinstance(const DT *C, const IT &n, const IT &m, DT *Cext, DT valmx = std::numeric_limits<DT>::max())
{
  IT msub = m-1, nsub = n-1, i, j;
  IT mpn = msub+nsub;
  // copy subsitutions
  for (j = 0; j < msub; j++)
    std::memcpy(Cext+j*(mpn),C+j*n,sizeof(DT)*nsub);
  // copy insertions
  for (i = nsub; i < mpn; i++)
    for (j = 0; j < msub; j++)
      if (i != j) Cext[j*mpn+i] = valmx;
      else Cext[j*mpn+i] = C[j*n+nsub];
  // copy removals
  for (i = 0; i < nsub; i++)
    for (j = msub; j < mpn; j++)
      if (i != j) Cext[j*mpn+i] = valmx;
      else Cext[j*mpn+i] = C[msub*n+i];
  // set completness
  for (i = nsub; i < mpn; i++)
    for (j = msub; j < mpn; j++)
      Cext[j*mpn+i] = 0;
}

// -----------------------------------------------------------
#endif
