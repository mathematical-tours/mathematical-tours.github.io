// =========================================================================
/** \file hungarian-lsape.hh
 *  \brief Hungarian algorithm for solving the Linear Sum Assignment Problem with Error-correction (LSAPE),
 *  aka minimal-cost error-correcting bipartite graph matching,
 * and its dual problem, according to a given edit cost matrix
 * \author Sebastien Bougleux (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072)
 * \author Luc Brun (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072)
 */
// =========================================================================
/* Hungarian algorithm for solving the Linear Sum Assignment
   Problem with Edition (LSAPE)
   
   authors: Sebastien Bougleux and Luc Brun
   institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072 
   
   -----------------------------------------------------------
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.

   -----------------------------------------------------------
   content
   -----------------------------------------------------------
   void hungarianLSAPE(C,n+1,m+1,rho,varrho,u,v,forb_assign)

   Compute a solution to the LSAPE.
   Compute an assignment with error-correction between two sets
   U={0,...,n-1} and V={0,...,m-1}, provided 
   a (n+1)x(m+1) edit cost matrix C, and such that the
   total cost sum of the assignment is miniminal. The 
   last row and the last column of C represent the costs
   of inserting and removing elements, respectively.

   we always have:
      C[n+(m+1)*m] = 0,  u[n] = v[m] = 0

   The resulting assignment is provided as two mappings
   rho:U->V\cup{m} and varrho:V->U\cup{n}.
   
   - Worst-case time complexity in O(min{n,m}^2max{n,m})
   - Space complexity in O(nm)

   reference: S. Bougleux and L. Brun
              Linear Sum Assignment with Edition
              Technical Report, March 2016
              Normandie Univ, GREYC UMR 6072
*/
// =========================================================================

#ifndef _HUNGARIAN_LSAPE_
#define _HUNGARIAN_LSAPE_

#include <limits>
#include <iostream>
// -----------------------------------------------------------
// Compute a initial partial assignment (rho,varrho) and associated dual variables (u,v)
// according to min on rows and then to min on reduced columns
// nass and mass are the number assigned elements in U and V respectively
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocLSAPE(const DT *C, const IT &nrows, const IT &ncols,
		       IT *rho, IT *varrho, DT *u, DT *v, IT &nass, IT &mass)
{
  const IT n = nrows-1, m = ncols-1;
  IT i = 0, j;
  DT mn, val;
  nass = mass = 0;
  u[n] = v[m] = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < ncols; j++)
    {
      const DT &c = C[j*nrows+i];
      if (c < mn) mn = c;
    }
    rho[i] = -1;
    u[i] = mn;
  }
  
  // find de min of each column
  for (j = 0; j < m; j++)
  {
    mn = std::numeric_limits<DT>::max();
    for (i = 0; i < nrows; i++)
    {
      val = C[j*nrows+i] - u[i];
      if (val < mn) mn = val;
    }
    v[j] = mn;
    varrho[j] = -1;
  }
  
  // assign
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      if (rho[i] < 0 && varrho[j] < 0 && C[j*nrows+i] == u[i] + v[j])
	{ rho[i] = j; varrho[j] = i; nass++; mass++; break; }
    }
    // last column
    if (rho[i] < 0 && C[m*nrows+i] == u[i]) { rho[i] = m; nass++; }
  }
  
  // last row
  for (j = 0; j < m; j++)
    if (varrho[j] < 0 && C[j*nrows+n] == v[j]) { varrho[j] = n; mass++; }
}

// -----------------------------------------------------------
// Compute a initial partial assignment (rho,varrho) and associated dual variables (u,v)
// according to min on rows and then to min on reduced columns
// nass and mass are the number assigned elements in U and V respectively
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocRowsLSAPE(const DT *C, const IT &nrows, const IT &ncols,
			   IT *rho, IT *varrho, DT *u, DT *v, IT &nass, IT &mass)
{
  const IT n = nrows-1, m = ncols-1;
  IT i = 0, j;
  DT mn, val;
  nass = mass = 0;
  u[n] = v[m] = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < ncols; j++)
    {
      const DT &c = C[j*nrows+i];
      if (c < mn) mn = c;
    }
    rho[i] = -1;
    u[i] = mn;
  }
  
  // find de min of each column
  for (j = 0; j < m; j++) { v[j] = 0; varrho[j] = -1; }
  
  // assign
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      if (rho[i] < 0 && varrho[j] < 0 && C[j*nrows+i] == u[i] + v[j])
	{ rho[i] = j; varrho[j] = i; nass++; mass++; break; }
    }
    // last column
    if (rho[i] < 0 && C[m*nrows+i] == u[i]) { rho[i] = m; nass++; }
  }
  
  // last row
  for (j = 0; j < m; j++)
    if (varrho[j] < 0 && C[j*nrows+n] == v[j]) { varrho[j] = n; mass++; }
}

// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocColsLSAPE(const DT *C, const IT &nrows, const IT &ncols,
			   IT *rho, IT *varrho, DT *u, DT *v, IT &nass, IT &mass)
{
  const IT n = nrows-1, m = ncols-1;
  IT i = 0, j = 0;
  DT mn, val;
  nass = mass = 0;
  u[n] = v[m] = 0;
  
  // find the min of each col
  for (; j < m; j++)
  {
    mn = std::numeric_limits<DT>::max();
    for (i = 0; i < nrows; i++)
    {
      const DT &c = C[j*nrows+i];
      if (c < mn) mn = c;
    }
    varrho[j] = -1;
    v[j] = mn;
  }
  
  // find de min of each row
  for (i = 0; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < ncols; j++)
    {
      val = C[j*nrows+i] - v[j];
      if (val < mn) mn = val;
    }
    u[i] = mn;
    rho[i] = -1;
  }
  
  // assign
  for (j = 0; j < m; j++)
  {
    for (i = 0; i < n; i++)
    {
      if (rho[i] < 0 && varrho[j] < 0 && C[j*nrows+i] == u[i] + v[j])
	{ rho[i] = j; varrho[j] = i; nass++; mass++; break; }
    }
    // last row n
    if (varrho[j] < 0 && C[j*nrows+n] == v[j]) { varrho[j] = n; mass++; }
  }
  
  // last column
  for (i = 0; i < n; i++)
    if (rho[i] < 0 && C[m*nrows+i] == u[i]) { rho[i] = m; nass++; }
}

// -----------------------------------------------------------
// same as above but with forbidden assignments
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocLSAPE_forb(const DT *C, const IT &nrows, const IT &ncols,
			   IT *rho, IT *varrho, DT *u, DT *v, IT &nass, IT &mass)
{
  const IT n = nrows-1, m = ncols-1;
  IT i = 0, j;
  DT mn, val;
  nass = mass = 0;
  u[n] = v[m] = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < ncols; j++)
    {
      const DT &c = C[j*nrows+i];
      if (c < 0) continue;
      if (c < mn) mn = c;
    }
    rho[i] = -1;
    u[i] = mn;
  }
  
  // find de min of each column
  for (j = 0; j < m; j++)
  {
    mn = std::numeric_limits<DT>::max();
    for (i = 0; i < nrows; i++)
    {
      const DT &c = C[j*nrows+i];
      if (c < 0) continue;
      val = c - u[i];
      if (val < mn) mn = val;
    }
    v[j] = mn;
    varrho[j] = -1;
  }
  
  // assign
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      const DT &c = C[j*nrows+i];
      if (c < 0) continue;
      if (rho[i] < 0 && varrho[j] < 0 && c == u[i] + v[j])
	{ rho[i] = j; varrho[j] = i; nass++; mass++; break; }
    }
    // last column
    if (rho[i] < 0 && C[m*nrows+i] == u[i]) { rho[i] = m; nass++; }
  }
  
  // last row
  for (j = 0; j < m; j++)
    if (varrho[j] < 0 && C[j*nrows+n] == v[j]) { varrho[j] = n; mass++; }
}

// -----------------------------------------------------------
// Construct an alternating tree rooted in an element of V
// and perform dual updates until an augmenting path is found
// -----------------------------------------------------------
template<class DT, typename IT>
void augmentingPathColLSAPE(const IT &k, const DT *C, const IT &nrows, const IT &m, const IT *rho, const IT *varrho, 
			    IT *U, IT *SV, IT *pred, DT *u, DT *v, DT *pi, IT &zi, IT &zj)
{
  const IT n = nrows-1;
  IT i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
  DT delta = 0, cred = 0;
  const IT *svptr = NULL, *luptr = NULL, *uend = U+n, *lusutop = U;
  bool lstrw = false;
  *SV = -1;
  zj = zi = -1;
  
  for (i = 0; i < n; i++) { pi[i] = std::numeric_limits<DT>::max(); U[i] = i; }
  
  while (true)
  {
    *SVptr = j; *(++SVptr) = -1;
    // last row: null element epsilon, it is a sink of an alternating path
    if (varrho[j] < n && C[j*nrows+n] == v[j]) { zi = n; zj = j; return; }
    
    for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
    {
      r = *uluptr;
      cred = C[j*nrows+r] - (u[r] + v[j]);
      if (cred < pi[r])
      {
	pred[r] = j;
	pi[r] = cred;
	if (cred == 0)
	{
	  if (rho[r] == -1 || rho[r] == m) { zi = r; zj = -1; return; }
	  i = *ulutop; *ulutop = r; *uluptr = i; ++ulutop;
	}
      }
    }
    
    if (lusutop == ulutop) // dual update
    {
      delta = std::numeric_limits<DT>::max(); lstrw = false;
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	if (pi[*uluptr] < delta) delta = pi[*uluptr];
      for (svptr = SV; *svptr != -1; ++svptr) // last row
      {
	cred = C[*svptr*nrows+n] - v[*svptr];
	if (cred <= delta) { delta = cred; lstrw = true; zj = *svptr; }
      }
      for (svptr = SV; *svptr != -1; ++svptr) v[*svptr] += delta;
      for (luptr = U; luptr != ulutop; ++luptr) u[*luptr] -= delta;
      if (lstrw) { zi = n; return; }
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
      {
	pi[*uluptr] -= delta;
	if (pi[*uluptr] == 0) 
	{
	  if (rho[*uluptr] == -1 || rho[*uluptr] == m) { zi = *uluptr; zj = -1; return; }
	  r = *ulutop; *ulutop = *uluptr; *uluptr = r; ++ulutop;
	}
      }
    } // end dual update
    i = *lusutop; ++lusutop; // i is now in SU
    j = rho[i];
  }
}

// -----------------------------------------------------------
// Construct an alternating tree rooted in an element of U
// and perform dual updates until an augmenting path is found
// -----------------------------------------------------------
template<class DT, typename IT>
void augmentingPathRowLSAPE(const IT &k, const DT *C, const IT &nrows, const IT &m, const IT *rho, const IT *varrho, 
			    IT *V, IT *SU, IT *pred, DT *u, DT *v, DT *pi, IT &zi, IT &zj)
{
  const IT n = nrows-1, *suptr = NULL, *lvptr = NULL, *vend = V+m, *lvsvtop = V;
  IT i = k, j = 0, c = 0, *SUptr = SU, *vlvtop = V, *vlvptr = NULL;
  DT delta = 0, cred = 0;
  bool lstcl = false;
  *SU = -1;
  zj = zi = -1;
  
  for (j = 0; j < m; j++) { pi[j] = std::numeric_limits<DT>::max(); V[j] = j; }
  
  while (true)
  {
    *SUptr = i; *(++SUptr) = -1;
    
    // last column: null element epsilon, it is a sink of an alternating path
    if (rho[i] < m && C[m*nrows+i] == u[i]) { zi = i; zj = m; return; }
    
    for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // U\LU
    {
      c = *vlvptr;
      cred = C[c*nrows+i] - (u[i] + v[c]);
      if (cred < pi[c])
      {
	pred[c] = i;
	pi[c] = cred;
	if (cred == 0)
	{
	  if (varrho[c] == -1 || varrho[c] == n) { zi = -1; zj = c; return; }
	  j = *vlvtop; *vlvtop = c; *vlvptr = j; ++vlvtop;
	}
      }
    }
        
    if (lvsvtop == vlvtop) // dual update
    {
      delta = std::numeric_limits<DT>::max(); lstcl = false;
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	if (pi[*vlvptr] < delta) delta = pi[*vlvptr];
      for (suptr = SU; *suptr != -1; ++suptr) // last column
      {
	cred = C[m*nrows+*suptr] - u[*suptr];
	if (cred <= delta) { delta = cred; lstcl = true; zi = *suptr; }
      }
      for (suptr = SU; *suptr != -1; ++suptr) u[*suptr] += delta;
      for (lvptr = V; lvptr != vlvtop; ++lvptr) v[*lvptr] -= delta;
      if (lstcl) { zj = m; return; }
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	{
	  pi[*vlvptr] -= delta;
	  if (pi[*vlvptr] == 0)
	  {
	    if (varrho[*vlvptr] == -1 || varrho[*vlvptr] == n) { zi = -1; zj = *vlvptr; return; }
	    c = *vlvtop; *vlvtop = *vlvptr; *vlvptr = c; ++vlvtop;
	  }
	}
    } // end dual update
    j = *lvsvtop; ++lvsvtop; // j is now in SV
    //if (varrho[j] == -1 || varrho[j] == n) { zi = -1; zj = j; break; }
    i = varrho[j];
  }
}

// -----------------------------------------------------------
// Construct an alternating tree rooted in an element of V
// and perform dual updates until an augmenting path is found
// -----------------------------------------------------------
template<class DT, typename IT>
void augmentingPathColLSAPE_forb(const IT &k, const DT *C, const IT &nrows, const IT &m, const IT *rho, const IT *varrho, 
				 IT *U, IT *SV, IT *pred, DT *u, DT *v, DT *pi, IT &zi, IT &zj)
{
  const IT n = nrows-1;
  IT i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
  DT delta = 0, cred = 0;
  const IT *svptr = NULL, *luptr = NULL, *uend = U+n, *lusutop = U;
  bool lstrw = false;
  *SV = -1;
  zj = zi = -1;
  for (i = 0; i < n; i++) { pi[i] = std::numeric_limits<DT>::max(); U[i] = i; }
  
  while (true)
  {
    *SVptr = j; *(++SVptr) = -1;
    
    // last row: null element epsilon, it is a sink of an alternating path
    if (varrho[j] < n && C[j*nrows+n] == v[j]) { zi = n; zj = j; return; }
    
    for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
    {
      r = *uluptr;
      if (C[j*nrows+r] < 0) continue; // check constraints on C (forbidden assignments)
      cred = C[j*nrows+r] - (u[r] + v[j]);
      if (cred < pi[r])
      {
	pred[r] = j;
	pi[r] = cred;
	if (cred == 0)
	{
	  if (rho[r] == -1 || rho[r] == m) { zi = r; zj = -1; return; }
	  i = *ulutop; *ulutop = r; *uluptr = i; ++ulutop;
	}
      }
    }
    
    if (lusutop == ulutop) // dual update
    {
      delta = std::numeric_limits<DT>::max(); lstrw = false;
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	if (pi[*uluptr] < delta) delta = pi[*uluptr];
      for (svptr = SV; *svptr != -1; ++svptr) // last row
      {
	cred = C[*svptr*nrows+n] - v[*svptr];
	if (cred <= delta) { delta = cred; lstrw = true; zj = *svptr; }
      }
      for (svptr = SV; *svptr != -1; ++svptr) v[*svptr] += delta;
      for (luptr = U; luptr != ulutop; ++luptr) u[*luptr] -= delta;
      if (lstrw) { zi = n; return; }
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
      {
	pi[*uluptr] -= delta;
	if (pi[*uluptr] == 0)
	{
	  if (rho[*uluptr] == -1 || rho[*uluptr] == m) { zi = *uluptr; zj = -1; return; }
	  r = *ulutop; *ulutop = *uluptr; *uluptr = r; ++ulutop;
	}
      }
    } // end dual update
    i = *lusutop; ++lusutop; // i is now in SU
    j = rho[i];
  }
}

// -----------------------------------------------------------
// Construct an alternating tree rooted in an element of U
// and perform dual updates until an augmenting path is found
// -----------------------------------------------------------
template<class DT, typename IT>
void augmentingPathRowLSAPE_forb(const IT &k, const DT *C, const IT &nrows, const IT &m, const IT *rho, const IT *varrho, 
				 IT *V, IT *SU, IT *pred, DT *u, DT *v, DT *pi, IT &zi, IT &zj)
{
  const IT n = nrows-1, *suptr = NULL, *lvptr = NULL, *vend = V+m, *lvsvtop = V;
  IT i = k, j = 0, c = 0, *SUptr = SU, *vlvtop = V, *vlvptr = NULL;
  DT delta = 0, cred = 0;
  bool lstcl = false;
  *SU = -1;
  zj = zi = -1;
  
  for (j = 0; j < m; j++) { pi[j] = std::numeric_limits<DT>::max(); V[j] = j; }
  
  while (true)
  {
    *SUptr = i; *(++SUptr) = -1;
    
    // last column: null element epsilon, it is a sink of an alternating path
    if (rho[i] < m && C[m*nrows+i] == u[i]) { zi = i; zj = m; return; }
    
    for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // U\LU
    {
      c = *vlvptr;
      if (C[c*nrows+i] < 0) continue; // check constraints on C (forbidden assignments)
      cred = C[c*nrows+i] - (u[i] + v[c]);
      if (cred < pi[c])
      {
	pred[c] = i;
	pi[c] = cred;
	if (cred == 0)
	{
	  if (varrho[c] == -1 || varrho[c] == n) { zi = -1; zj = c; return; }
	  j = *vlvtop; *vlvtop = c; *vlvptr = j; ++vlvtop;
	}
      }
    }
        
    if (lvsvtop == vlvtop) // dual update
    {
      delta = std::numeric_limits<DT>::max(); lstcl = false;
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	if (pi[*vlvptr] < delta) delta = pi[*vlvptr];
      for (suptr = SU; *suptr != -1; ++suptr) // last column
      {
	cred = C[m*nrows+*suptr] - u[*suptr];
	if (cred <= delta) { delta = cred; lstcl = true; zi = *suptr; }
      }
      for (suptr = SU; *suptr != -1; ++suptr) u[*suptr] += delta;
      for (lvptr = V; lvptr != vlvtop; ++lvptr) v[*lvptr] -= delta;
      if (lstcl) { zj = m; return; }
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	{
	  pi[*vlvptr] -= delta;
	  if (pi[*vlvptr] == 0)
	  {
	    if (varrho[*vlvptr] == -1 || varrho[*vlvptr] == n) { zi = -1; zj = *vlvptr; return; }
	    c = *vlvtop; *vlvtop = *vlvptr; *vlvptr = c; ++vlvtop;
	  }
	}
    } // end dual update
    j = *lvsvtop; ++lvsvtop; // j is now in SV
    i = varrho[j];
  }
}

// -----------------------------------------------------------
// Main function: Hungarian algorithm for LSAPE
// -----------------------------------------------------------
/**
 * \brief Compute a solution to the LSAPE (minimal-cost error-correcting bipartite graph matching) with the Hungarian method
 * \param[in] C nrowsxncols edit cost matrix represented as an array if size \p nrows.ncols obtained by concatenating its columns, column \p nrows-1 are the costs of removing the elements of the 1st set, and the row \p ncols-1 represents the costs of inserting an element of the 2nd set
 * \param[in] nrows Number of rows of \p C
 * \param[in] ncols Number of columns of \p C
 * \param[out] rho Array of size \p nrows-1 (must be previously allocated), rho[i]=j indicates that i is assigned to j (substituted by j if j<ncols-1, or removed if j=ncols-1)
 * \param[out] varrho Array of size \p m (must be previously allocated), varrho[j]=i indicates that j is assigned to i (substituted to i if i<nrows-1, or inserted if i=nrows)
 * \param[out] u Array of dual variables associated to the 1st set (rows of \p C), of size \p nrows
 * \param[out] v Array of dual variables associated to the 2nd set (columns of \p C), of size \p ncols
 * \param[in] forb_assign If true, forbidden assignments are marked with negative values in the cost matrix
 * \details A solution to the LSAPE is computed with the primal-dual version of the Hungarian algorithm, as detailed in:
 * \li <em>S. Bougleux and L. Brun, Linear Sum Assignment with Edition, Technical Report, Normandie Univ, GREYC UMR 6072, 2016</em>
 *
 * This version updates dual variables \c u and \c v, and at each iteration, the current matching is augmented by growing only one Hungarian tree until an augmenting path is found. Our implementation uses a Bread-First-like strategy to construct the tree, according to a FIFO strategy to select the next element at each iteration of the growing process.
 *
 * Complexities:
 * \li O(min{n,m}²max{n,m}) in time (worst-case)
 * \li O(nm) in space
 *
 * \remark
 * Template \p DT allows to compute a solution with integer or floating-point values. Note that rounding errors may occur with floating point values when dual variables are updated but this does not affect the overall process.
 */
template <class DT, typename IT>
void hungarianLSAPE(const DT *C, const IT &nrows, const IT &ncols, IT *rho, IT *varrho, DT *u, DT *v, unsigned short init_type = 1, bool forb_assign = false)
{
  const IT n = nrows-1, m = ncols-1;
  IT nass = 0, mass = 0, i = -1, j = -1, r = -1, c = -1, k = -1, *S = NULL, *U = NULL, *V = NULL, *pred = NULL;
  DT *pi = NULL;
  u[n] = v[m] = 0;

  if (init_type == 0)
  { for (k = 0; k < n; k++) { rho[k] = -1; u[k] = 0; } for (k = 0; k < m; k++) { varrho[k] = -1; v[k] = 0; } }
  u[n] = v[m] = 0;
  
  if (forb_assign)
  {
    // initialization -----------------------------------------------
    if (init_type == 1) basicPreprocLSAPE_forb<DT,IT>(C,nrows,ncols,rho,varrho,u,v,nass,mass);
    
    // augmentation of columns --------------------------------------
    if (mass < m)
    {
      U = new IT[nrows]; S = new IT[ncols];  pi = new DT[n]; pred = new IT[n];
      for (k = 0; k < m; k++)
	if (varrho[k] == -1)
	{
	  augmentingPathColLSAPE_forb<DT,IT>(k,C,nrows,m,rho,varrho,U,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
	  if (i == n) { r = varrho[j]; varrho[j] = i; i = r; }
	  else j = -1;
	  for (; j != k;)  // update primal solution = new partial assignment
	  {
	    j = pred[i]; rho[i] = j;
	    r = varrho[j]; varrho[j] = i; i = r;
	  }
	}
      delete[] U; delete[] pred; delete[] pi; delete[] S;
    }
    // augmentation of rows --------------------------------------
    if (nass < n)
    {
      V = new IT[ncols]; S = new IT[nrows];  pi = new DT[m]; pred = new IT[m];
      for (k = 0; k < n; k++)
	if (rho[k] == -1)
	{
	  augmentingPathRowLSAPE_forb<DT,IT>(k,C,nrows,m,rho,varrho,V,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
	  if (j == m) { c = rho[i]; rho[i] = j; j = c; }
	  else i = -1;
	  for (; i != k;)  // update primal solution = new partial assignment
	  {
	    i = pred[j]; varrho[j] = i;
	    c = rho[i]; rho[i] = j; j = c;
	  }
	}
      delete[] V; delete[] pred; delete[] pi; delete[] S;
    }
  }
  else
  {
    //if (nrows < ncols)
    //{
      // initialization -----------------------------------------------
      if (init_type == 1) basicPreprocColsLSAPE<DT,IT>(C,nrows,ncols,rho,varrho,u,v,nass,mass);
      // augmentation of columns --------------------------------------
      if (mass < m)
      {
	U = new IT[nrows]; S = new IT[ncols];  pi = new DT[n]; pred = new IT[n];
	for (k = 0; k < m; k++)
	  if (varrho[k] == -1)
	  {
	    augmentingPathColLSAPE<DT,IT>(k,C,nrows,m,rho,varrho,U,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
	    if (i == n) { r = varrho[j]; varrho[j] = i; i = r; }
	    else j = -1;
	    for (; j != k;)  // update primal solution = new partial assignment
	    {
	      j = pred[i]; rho[i] = j;
	      r = varrho[j]; varrho[j] = i; i = r;
	    }
	  }
	delete[] U; delete[] pred; delete[] pi; delete[] S;
      }
      // augmentation of rows --------------------------------------
      if (nass < n)
      {
	V = new IT[ncols]; S = new IT[nrows];  pi = new DT[m]; pred = new IT[m];
	for (k = 0; k < n; k++)
	  if (rho[k] == -1)
	  {
	    augmentingPathRowLSAPE<DT,IT>(k,C,nrows,m,rho,varrho,V,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
	    if (j == m) { c = rho[i]; rho[i] = j; j = c; }
	    else i = -1;
	    for (; i != k;)  // update primal solution = new partial assignment
	    {
	      i = pred[j]; varrho[j] = i;
	      c = rho[i]; rho[i] = j; j = c;
	    }
	  }
	delete[] V; delete[] pred; delete[] pi; delete[] S;
      }
      //}
    // else
    // {
    //   // augmentation of rows --------------------------------------
    //   if (nass < n)
    //   {
    // 	V = new IT[ncols]; S = new IT[nrows];  pi = new DT[m]; pred = new IT[m];
    // 	for (k = 0; k < n; k++)
    // 	  if (rho[k] == -1)
    // 	  {
    // 	    augmentingPathRowLSAPE<DT,IT>(k,C,nrows,m,rho,varrho,V,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
    // 	    if (j == m) { c = rho[i]; rho[i] = j; j = c; }
    // 	    else i = -1;
    // 	    for (; i != k;)  // update primal solution = new partial assignment
    // 	    {
    // 	      i = pred[j]; varrho[j] = i;
    // 	      c = rho[i]; rho[i] = j; j = c;
    // 	    }
    // 	  }
    // 	delete[] V; delete[] pred; delete[] pi; delete[] S;
    //   }
    //   // augmentation of columns --------------------------------------
    //   if (mass < m)
    //   {
	
    // 	U = new IT[nrows]; S = new IT[ncols];  pi = new DT[n]; pred = new IT[n];
    // 	for (k = 0; k < m; k++)
    // 	  if (varrho[k] == -1)
    // 	  {
    // 	    //	  std::cerr << "yo\n";
    // 	    augmentingPathColLSAPE<DT,IT>(k,C,nrows,m,rho,varrho,U,S,pred,u,v,pi,i,j);  // augment always finds an augmenting path
    // 	    if (i == n) { r = varrho[j]; varrho[j] = i; i = r; }
    // 	    else j = -1;
    // 	    for (; j != k;)  // update primal solution = new partial assignment
    // 	    {
    // 	      j = pred[i]; rho[i] = j;
    // 	      r = varrho[j]; varrho[j] = i; i = r;
    // 	    }
    // 	  }
    // 	delete[] U; delete[] pred; delete[] pi; delete[] S;
    //   }
      // }
  }
}

// -----------------------------------------------------------
// Reconstruct the list of insertions epsilon->V
// return false if no insertion is found
// -----------------------------------------------------------
template <typename IT>
IT reconstructInsertions(const IT *varrho, const IT &n, const IT &m, IT **rhoeps)
{
  IT nb = 0;
  for (IT j = 0; j < m; j++) if (varrho[j] == n) nb++;
  if (nb > 0) *rhoeps = new IT[nb];
  else return nb;
  for (IT j = 0, k = 0; j < m; j++) if (varrho[j] == n) { (*rhoeps)[k] = j; k++; }
  return nb;
}

// -----------------------------------------------------------
#endif
