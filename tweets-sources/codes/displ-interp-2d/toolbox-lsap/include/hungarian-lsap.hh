// =========================================================================
/** \file hungarian-lsap.hh
 *  \brief Hungarian algorithm for solving symmetric or assymetric Linear Sum Assignment Problem (LSAP), 
 *  and its dual problem
 * \author Sebastien Bougleux (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France)
 * \author Luc Brun (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France)
 */
// =========================================================================
/* 
   -----------------------------------------------------------
   
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.
   
   You should have received a copy of the CeCILL-C License along with 
   LSAPE. If not, see http://www.cecill.info/index.en.html

   -----------------------------------------------------------
   
   Creation: December 5 2015
   Modifications: July 5 2017
   
   -----------------------------------------------------------
   Content
   -----------------------------------------------------------
   
   void hungarianLSAP(C,n,m,rho,u,v,forbidden_assign=false)
   
   Compute a solution to the linear sum assignment problem (LSAP) between two sets U and V,
   as well as its dual problem (element labeling), provided a nxm cost matrix C
   
   Inputs:
   - nxm matrix C of non-negative costs
     -> C is assumed to be already allocated as an array of size nm, obtained by concatenating columns (as in matlab)
   - optional boolean forbidden_assign (false by default)
     -> if true, indicates that forbidden assignments are marked with negative values in the cost matrix
     -> in this case, these assignments are not treated
     -> if false, the matrix C should not contain negative values
     
     for example, the following matrix contains forbidden assignments:
          _          _
         |  10  5  4  |
     C = |  -1  2  7  |, it is represented as the array [10 -1 5 5 2 -1 4 7 8]
         |_  5 -1  8 _|
         
   Outputs:
   - when n=m, the resulting optimal assignment/perfect matching is provided as a permutation of {0,...,n-1}, represented by the array rho of size n
   - else, the resulting optimal assignment/maximum matching is provided as a partial permutation of {0,...,m-1}, represented by the array rho of size n
   - (u,v) is the dual optimal labeling (the labels of the elements in each set), each dual variable is represented by an array of size n
     -> u[i] is the label of element i in the first set
     -> v[j] is the label of element j in the 2nd set

   Remarks:
   - the minimal cost can be computed by summing the values of the optimal labeling (u,v): sum { u[i]+v[i], i=1,...,n }
     -> you can use the function labelingSum(u,v) to obtain this result
   
   Complexity (worst-case):
   - 0(max{n,m}min{n,m}²) in time
   - 0(nm) in space

   References:
   - E.L. Lawler. Combinatorial Optimization: Networks and Matroids. Holt, Rinehart and Winston, New York, 1976
   - R. Burkard, M. Dell’Amico and S. Martello. Assignment Problems. SIAM, 2009
   
   Remark: a FIFO strategy is used to construct Hungarian trees and so to find augmenting paths
*/
// =========================================================================

#ifndef _HUNGARIAN_LSAP_
#define _HUNGARIAN_LSAP_

#include <limits>

// -----------------------------------------------------------
// Basic preprocessing for rows
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocRectRowsLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho, DT *u, DT *v, IT &nass)
{
  IT i = 0, j;
  DT mn;
  nass = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < m; j++)
    {
      const DT &c = C[j*n+i];
      if (c < mn) mn = c;
    }
    u[i] = mn;
    rho[i] = -1;
  }
  
  for (j = 0; j < m; j++) { v[j] = 0; varrho[j] = -1; }
  
  // assign
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      if (rho[i] == -1 && varrho[j] == -1 && C[j*n+i] == u[i])
	{ rho[i] = j; varrho[j] = i; nass++; break; }
    }
  }
}

// -----------------------------------------------------------
// Basic preprocessing for columns
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocRectColsLSAP(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho, DT *u, DT *v, IT &mass)
{
  IT i, j = 0;
  DT mn;
  mass = 0;
  
  // find the min of each row
  for (; j < m; j++)
  {
    mn = std::numeric_limits<DT>::max();
    for (i = 0; i < n; i++)
    {
      const DT &c = C[j*n+i];
      if (c < mn) mn = c;
    }
    v[j] = mn;
    varrho[j] = -1;
  }
  
  for (i = 0; i < n; i++) { u[i] = 0; rho[i] = -1; }
  
  // assign
  for (j = 0; j < m; j++)
  {
    for (i = 0; i < n; i++)
    {
      if (rho[i] == -1 && varrho[j] == -1 && C[j*n+i] == v[j])
	{ rho[i] = j; varrho[j] = i; mass++; break; }
    }
  }
}

// -----------------------------------------------------------
// Basic preprocessing with forbidden assignments for rows
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocRectRowsLSAP_forb(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho, DT *u, DT *v, IT &nass)
{
  IT i = 0, j;
  DT mn;
  nass = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = std::numeric_limits<DT>::max();
    for (j = 0; j < m; j++)
    {
      const DT &c = C[j*n+i];
      if (c < 0) continue;
      if (c < mn) mn = c;
    }
    u[i] = mn;
    rho[i] = -1;
  }
  
  for (j = 0; j < m; j++) { v[j] = 0; varrho[j] = -1; }
  
  // assign
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      const DT &c = C[j*n+i];
      if (c < 0) continue;
      if (rho[i] == -1 && varrho[j] == -1 && c == u[i])
	{ rho[i] = j; varrho[j] = i; nass++; break; }
    }
  }
}

// -----------------------------------------------------------
// Basic preprocessing with forbidden assignments for columns
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocRectColsLSAP_forb(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho, DT *u, DT *v, IT &mass)
{
  IT i, j = 0;
  DT mn;
  mass = 0;
  
  // find the min of each row
  for (; j < m; j++)
  {
    mn = std::numeric_limits<DT>::max();
    for (i = 0; i < n; i++)
    {
      const DT &c = C[j*n+i];
      if (c < 0) continue;
      if (c < mn) mn = c;
    }
    v[j] = mn;
    varrho[j] = -1;
  }
  
  for (i = 0; i < n; i++) { u[i] = 0; rho[i] = -1; }
  
  // assign
  for (j = 0; j < m; j++)
  {
    for (i = 0; i < n; i++)
    {
      const DT &c = C[j*n+i];
      if (c < 0) continue;
      if (rho[i] == -1 && varrho[j] == -1 && c == v[j])
	{ rho[i] = j; varrho[j] = i; mass++; break; }
    }
  }
}

// -----------------------------------------------------------
// Construct an alternating (Hungarian) tree rooted in an uncovered element of V
// and stop when an uncovered element of U is reached
// -----------------------------------------------------------
template<class DT, typename IT>
void augmentingPathCol(const IT &k, const DT *C, const IT &n, const IT &m, const IT *rho, 
		       IT *U, IT *SV, IT *pred, DT *u, DT *v, DT *pi, IT &zi)
{
  IT i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
  DT delta = 0, mx = std::numeric_limits<DT>::max(), cred = 0, zero = 0;
  const IT *svptr = NULL, *luptr = NULL, *uend = U+n, *lusutop = U;
  *SV = -1;
  zi = -1;
  for (i = 0; i < n; i++) { pi[i] = mx; U[i] = i; }
  
  while (true)
  {
    *SVptr = j; *(++SVptr) = -1;
    for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
    {
      r = *uluptr;
      cred = C[j*n+r] - (u[r] + v[j]);
      if (cred < pi[r])
      {
	pred[r] = j;
	pi[r] = cred;
	if (cred == 0)
	{
	  if (rho[r] == -1) { zi = r; return; }
	  i = *ulutop; *ulutop = r; *uluptr = i; ++ulutop;
	}
      }
    }
    
    if (lusutop == ulutop) // dual update
    {
      //nb_dual_col++;
      delta = mx;
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	if (pi[*uluptr] < delta) delta = pi[*uluptr];
      for (svptr = SV; *svptr != -1; ++svptr) v[*svptr] += delta;
      for (luptr = U; luptr != ulutop; ++luptr) u[*luptr] -= delta;
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
      {
	pi[*uluptr] -= delta;
	if (pi[*uluptr] == 0) 
	{
	  if (rho[*uluptr] == -1) { zi = *uluptr; return; }
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
void augmentingPathRow(const IT &k, const DT *C, const IT &n, const IT &m, const IT *varrho, 
		       IT *V, IT *SU, IT *pred, DT *u, DT *v, DT *pi, IT &zj)
{
  const IT *suptr = NULL, *lvptr = NULL, *vend = V+m, *lvsvtop = V;
  IT i = k, j = 0, c = 0, *SUptr = SU, *vlvtop = V, *vlvptr = NULL;
  DT delta = 0, mx = std::numeric_limits<DT>::max(), cred = 0;
  *SU = -1;
  zj = -1;
  
  for (j = 0; j < m; j++) { pi[j] = mx; V[j] = j; }
  
  while (true)
  {
    *SUptr = i; *(++SUptr) = -1;
    for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // U\LU
    {
      c = *vlvptr;
      cred =  C[c*n+i] - (u[i] + v[c]);
      if (cred < pi[c])
      {
	pred[c] = i;
	pi[c] = cred;
	if (cred == 0)
	{
	  if (varrho[c] == -1) { zj = c; return; }
	  j = *vlvtop; *vlvtop = c; *vlvptr = j; ++vlvtop;
	}
      }
    }
        
    if (lvsvtop == vlvtop) // dual update
    {
      //nb_dual_row++;
      delta = mx;
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	if (pi[*vlvptr] < delta) delta = pi[*vlvptr];
      for (suptr = SU; *suptr != -1; ++suptr) u[*suptr] += delta;
      for (lvptr = V; lvptr != vlvtop; ++lvptr) v[*lvptr] -= delta;
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	{
	  pi[*vlvptr] -= delta;
	  if (pi[*vlvptr] == 0)
	  {
	    if (varrho[*vlvptr] == -1) { zj = *vlvptr; return; }
	    c = *vlvtop; *vlvtop = *vlvptr; *vlvptr = c; ++vlvtop;
	  }
	}
    } // end dual update
    j = *lvsvtop; ++lvsvtop; // j is now in SV
    i = varrho[j];
  }
}

// -----------------------------------------------------------
// Construct an alternating (Hungarian) tree rooted in an uncovered element of V
// and stop when an uncovered element of U is reached
// -----------------------------------------------------------
template<class DT, typename IT>
void augmentingPathCol_forb(const IT &k, const DT *C, const IT &n, const IT &m, const IT *rho,  
			    IT *U, IT *SV, IT *pred, DT *u, DT *v, DT *pi, IT &zi)
{
  IT i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
  DT delta = 0, mx = (DT)std::numeric_limits<DT>::max()/2, cred = 0, zero = 0, deux = 2;
  const IT *svptr = NULL, *luptr = NULL, *uend = U+n, *lusutop = U;
  *SV = -1;
  zi = -1;
  
  for (i = 0; i < n; i++) { pi[i] = mx; U[i] = i; }
  
  while (true)
  {
    *SVptr = j; *(++SVptr) = -1;
    
    for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
    {
      r = *uluptr;
      cred = C[j*n+r];
      if (cred < zero) continue;
      cred -= u[r] + v[j];
      if (cred < pi[r])
      {
	pred[r] = j;
	pi[r] = cred;
	if (cred == zero)
	{
	  if (rho[r] == -1) { zi = r; return; }
	  i = *ulutop; *ulutop = r; *uluptr = i; ++ulutop;
	}
      }
    }
    
    if (lusutop == ulutop) // dual update
    {
      delta = mx;
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	if (pi[*uluptr] < delta) delta = pi[*uluptr];
      for (svptr = SV; *svptr != -1; ++svptr) v[*svptr] += delta;
      for (luptr = U; luptr != ulutop; ++luptr) u[*luptr] -= delta;
      for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
      {
	pi[*uluptr] -= delta;
	if (pi[*uluptr] == zero) 
	{
	  if (rho[*uluptr] == -1) { zi = *uluptr; return; }
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
void augmentingPathRow_forb(const IT &k, const DT *C, const IT &n, const IT &m, const IT *varrho, 
			    IT *V, IT *SU, IT *pred, DT *u, DT *v, DT *pi, IT &zj)
{
  const IT *suptr = NULL, *lvptr = NULL, *vend = V+m, *lvsvtop = V;
  IT i = k, j = 0, c = 0, *SUptr = SU, *vlvtop = V, *vlvptr = NULL;
  DT delta = 0, mx = std::numeric_limits<DT>::max(), cred = 0;
  *SU = -1;
  zj = -1;
  
  for (j = 0; j < m; j++) { pi[j] = mx; V[j] = j; }
  
  while (true)
  {
    *SUptr = i; *(++SUptr) = -1;
    
    for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // U\LU
    {
      c = *vlvptr;
      cred = C[c*n+i];
      if (cred < 0) continue;
      cred -=  u[i] + v[c];
      if (cred < pi[c])
      {
	pred[c] = i;
	pi[c] = cred;
	if (cred == 0)
	{
	  if (varrho[c] == -1) { zj = c; return; }
	  j = *vlvtop; *vlvtop = c; *vlvptr = j; ++vlvtop;
	}
      }
    }
        
    if (lvsvtop == vlvtop) // dual update
    {
      delta = mx;
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	if (pi[*vlvptr] < delta) delta = pi[*vlvptr];
      for (suptr = SU; *suptr != -1; ++suptr) u[*suptr] += delta;
      for (lvptr = V; lvptr != vlvtop; ++lvptr) v[*lvptr] -= delta;
      for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	{
	  pi[*vlvptr] -= delta;
	  if (pi[*vlvptr] == 0)
	  {
	    if (varrho[*vlvptr] == -1) { zj = *vlvptr; return; }
	    c = *vlvtop; *vlvtop = *vlvptr; *vlvptr = c; ++vlvtop;
	  }
	}
    } // end dual update
    j = *lvsvtop; ++lvsvtop; // j is now in SV
    i = varrho[j];
  }
}

// --------------------------------------------------------------------------
// Main function for Hungarian algorithm for LSAP with rectangular cost matrices
// --------------------------------------------------------------------------
template <class DT, typename IT>
void hungarianRectLSAP(const DT *C, const IT &nrows, const IT &ncols, IT *rho, IT *varrho, DT *u, DT *v, unsigned short init_type = 1, bool forbidden_assign = false)
{
  IT nass = 0, mass = 0, i = -1, j = -1, r = -1, c = -1, k = -1, *S = NULL, *U = NULL, *V = NULL, *pred = NULL;
  DT *pi = NULL, mxC = std::numeric_limits<DT>::min();
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[ncols]; }

  if (forbidden_assign) // with forbidden assignments
  {
    if (nrows < ncols)
    {
      // initialization -----------------------------------------------
      switch(init_type)
      {
      case 0: { 
	for (IT i=0; i < nrows; i++) { u[i] = 0; rho[i] = -1; }
	for (IT j=0; j < ncols; j++) { u[j] = 0; varrho[j] = -1; }
      } break;
      case 1: basicPreprocRectRowsLSAP_forb(C,nrows,ncols,rho,varrho,u,v,nass);
      }
      
      // augmentation of rows --------------------------------------
      if (nass < nrows)
      {
	V = new IT[ncols+1]; S = new IT[nrows+1];  pi = new DT[ncols]; pred = new IT[ncols];
	for (k = 0; k < nrows; k++)
	  if (rho[k] == -1)
	  {
	    augmentingPathRow_forb(k,C,nrows,ncols,varrho,V,S,pred,u,v,pi,j);  // augment always finds an augmenting path
	    // augment
	    for (i = -1; i != k;)  // update primal solution = new partial assignment
	    {
	      i = pred[j]; varrho[j] = i;
	      c = rho[i]; rho[i] = j; j = c;
	    }
	    nass++;
	  }
	delete[] V; delete[] pred; delete[] pi; delete[] S;
      }
    }
    else // nrows > ncols
    {
      // initialization -----------------------------------------------
      switch(init_type)
      {
      case 0: { 
	for (IT i=0; i < nrows; i++) { u[i] = 0; rho[i] = -1; }
	for (IT j=0; j < ncols; j++) { u[j] = 0; varrho[j] = -1; }
      } break;
      case 1: basicPreprocRectColsLSAP_forb(C,nrows,ncols,rho,varrho,u,v,mass);
      }
    
      // augmentation of columns --------------------------------------
      if (mass < ncols)
      {
	U = new IT[nrows+1]; S = new IT[ncols+1];  pi = new DT[nrows]; pred = new IT[nrows];
	for (k = 0; k < ncols; k++)
	  if (varrho[k] == -1)
	  {
	    augmentingPathCol_forb(k,C,nrows,ncols,rho,U,S,pred,u,v,pi,i);  // augment always finds an augmenting path
	    for (j = -1; j != k;)  // update primal solution = new partial assignment
	    {
	      j = pred[i]; rho[i] = j;
	      r = varrho[j]; varrho[j] = i; i = r;
	    }
	    mass++;
	  }
	delete[] U; delete[] pred; delete[] pi; delete[] S;
      }
    }
  }
  else // without forbidden assignments
  {
    if (nrows < ncols)
    {
      // initialization -----------------------------------------------
      switch(init_type)
      {
      case 0: { 
	for (IT i=0; i < nrows; i++) { u[i] = 0; rho[i] = -1; }
	for (IT j=0; j < ncols; j++) { u[j] = 0; varrho[j] = -1; }
      } break;
      case 1: basicPreprocRectRowsLSAP(C,nrows,ncols,rho,varrho,u,v,nass);
      }

      // augmentation of rows --------------------------------------
      if (nass < nrows)
      {
	V = new IT[ncols+1]; S = new IT[nrows+1];  pi = new DT[ncols]; pred = new IT[ncols];
	for (k = 0; k < nrows; k++)
	  if (rho[k] == -1)
	  {
	    augmentingPathRow(k,C,nrows,ncols,varrho,V,S,pred,u,v,pi,j);  // augment always finds an augmenting path
	    // augment
	    for (i = -1; i != k;)  // update primal solution = new partial assignment
	    {
	      i = pred[j]; varrho[j] = i;
	      c = rho[i]; rho[i] = j; j = c;
	    }
	    nass++;
	  }
	delete[] V; delete[] pred; delete[] pi; delete[] S;
      }
    }
    else // nrows > ncols
    {
      // initialization -----------------------------------------------
      switch(init_type)
      {
      case 0: { 
	for (IT i=0; i < nrows; i++) { u[i] = 0; rho[i] = -1; }
	for (IT j=0; j < ncols; j++) { u[j] = 0; varrho[j] = -1; }
      } break;
      case 1: basicPreprocRectColsLSAP(C,nrows,ncols,rho,varrho,u,v,mass);
      }
    
      // augmentation of columns --------------------------------------
      if (mass < ncols)
      {
	U = new IT[nrows+1]; S = new IT[ncols+1];  pi = new DT[nrows]; pred = new IT[nrows];
	for (k = 0; k < ncols; k++)
	  if (varrho[k] == -1)
	  {
	    augmentingPathCol(k,C,nrows,ncols,rho,U,S,pred,u,v,pi,i);  // augment always finds an augmenting path
	    for (j = -1; j != k;)  // update primal solution = new partial assignment
	    {
	      j = pred[i]; rho[i] = j;
	      r = varrho[j]; varrho[j] = i; i = r;
	    }
	    mass++;
	  }
	delete[] U; delete[] pred; delete[] pi; delete[] S;
      }
    }
  }
  if (deletevarrho) { delete[] varrho; varrho = NULL; }
}

// -----------------------------------------------------------
// Compute a initial partial assignment (rho,varrho) and associated dual variables (u,v)
// according to min on rows and then to min on reduced columns
// mass is the number elements assigned during the process
// -----------------------------------------------------------
template<class DT, typename IT>
void basicPreprocLSAP(const DT *C, const IT &n, IT *rho, IT *varrho, DT *u, DT *v, IT &mass)
{
  IT i = 0, j;
  DT mn, val, mx = std::numeric_limits<DT>::max();
  mass = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = mx;
    for (j = 0; j < n; j++)
    {
      const DT &c = C[j*n+i];
      if (c < mn) mn = c;
    }
    u[i] = mn;
    rho[i] = -1;
  }
  
  // find the min of each column
  for (j = 0; j < n; j++)
  {
    mn = mx;
    for (i = 0; i < n; i++)
    {
      val = C[j*n+i] - u[i];
      if (val < mn) mn = val;
    }
    v[j] = mn;
    varrho[j] = -1;
  }
  
  // assign
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (rho[i] == -1 && varrho[j] == -1 && C[j*n+i] == u[i] + v[j])
	{ rho[i] = j; varrho[j] = i; mass++; break; }
    }
  }
}

// -----------------------------------------------------------
// Compute a initial partial assignment (rho,varrho) and associated dual variables (u,v)
// version with negative values in C corresponding to forbidden assignments of elements
// according to min on rows and then to min on reduced columns
// mass is the number elements assigned during the process
// -----------------------------------------------------------
template <class DT, typename IT>
void basicPreprocLSAP_forb(const DT *C, const IT &n, IT *rho, IT *varrho, DT *u, DT *v, IT &mass)
{
  IT i = 0, j;
  DT mn, val, mx = std::numeric_limits<DT>::max();
  mass = 0;
  
  // find the min of each row
  for (; i < n; i++)
  {
    mn = mx;
    for (j = 0; j < n; j++)
    {
      const DT &c = C[j*n+i];
      if (c < 0) continue;
      if (c < mn) mn = c;
    }
    u[i] = mn;
    rho[i] = -1;
  }
  
  // find the min of each column
  for (j = 0; j < n; j++)
  {
    mn = mx;
    for (i = 0; i < n; i++)
    {
      const DT &c = C[j*n+i];
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
    for (j = 0; j < n; j++)
    {
      const DT &c = C[j*n+i];
      if (c < 0) continue;
      if (rho[i] == -1 && varrho[j] == -1 && c == u[i] + v[j])
	{ rho[i] = j; varrho[j] = i; mass++; break; }
    }
  }
}

// -----------------------------------------------------------
// Main Hungarian algorithm for LSAP with square cost matrices
// -----------------------------------------------------------
template <class DT, typename IT>
void hungarianSquareLSAP(const DT *C, const IT &n, IT *rho, DT *u, DT *v, IT *varrho = NULL, unsigned short init_type = 1, bool forbidden_assign = false)
{
  IT nass = 0, i = -1, j = -1, r = -1, k = -1;
  IT *SV = NULL, *pred = NULL, *U = NULL;
  DT *pi = NULL;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[n]; }

  if (forbidden_assign)
  {
    // initialization -----------------------------------------------
    switch (init_type)
    {
    case 0: { for (IT i=0; i < n; i++) { u[i] = v[i] = 0; rho[i] = varrho[i] = -1; } } break;
    case 1: basicPreprocLSAP_forb<DT,IT>(C,n,rho,varrho,u,v,nass);
    }

    // augmentation of columns --------------------------------------
    if (nass < n)
    {
      U = new IT[n+1]; SV = new IT[n+1];  pi = new DT[n+1]; pred = new IT[n];
      for (k = 0; k < n; k++)
	if (varrho[k] == -1)
	{
	  augmentingPathCol_forb<DT,IT>(k,C,n,n,rho,U,SV,pred,u,v,pi,i);  // augment always finds an augmenting path
	  for (j = -1; j != k;)  // update primal solution = new partial assignment
	  {
	    j = pred[i]; rho[i] = j;
	    r = varrho[j]; varrho[j] = i; i = r;
	  }
	}
      delete[] U; delete[] pred; delete[] pi; delete[] SV;
    }
  }
  // without forbidden assignments (classical case) -------------------------
  else
  {
    // initialization -----------------------------------------------
    switch (init_type)
    {
    case 0: { for (IT i=0; i < n; i++) { u[i] = v[i] = 0; rho[i] = varrho[i] = -1; } } break;
    case 1: basicPreprocLSAP<DT,IT>(C,n,rho,varrho,u,v,nass);
    }

    // augmentation of columns --------------------------------------
    if (nass < n)
    {
      U = new IT[n+1]; SV = new IT[n+1];  pi = new DT[n+1]; pred = new IT[n];
      for (k = 0; k < n; k++)
	if (varrho[k] == -1)
	{
	  augmentingPathCol<DT,IT>(k,C,n,n,rho,U,SV,pred,u,v,pi,i);  // augment always finds an augmenting path
	  for (j = -1; j != k;)  // update primal solution = new partial assignment
	  {
	    j = pred[i]; rho[i] = j;
	    r = varrho[j]; varrho[j] = i; i = r;
	  }
	}
      delete[] U; delete[] pred; delete[] pi; delete[] SV;
    }
  }
  if (deletevarrho) delete[] varrho;
}

// --------------------------------------------------------------------------------
// Main function: hungarian algorithm for both square and rectangular cost matrices
/**
 * \brief Compute a solution to symmetric or assymetric LSAP with the Hungarian algorithm
 * \param[in] C nxm cost matrix represented as an array of size \c nm obtained by concatenating its columns
 * \param[in] n Number of rows of \p C (size of the 1st set)
 * \param[in] m Number of columns of \p C (size of the 2nd set)
 * \param[out] rho A solution to the LSAP: an array of size \p n (must be previously allocated), rho[i]=-1 indicates that i is not assigned, else rho[i]=j indicates that i is assigned to j
 * \param[out] u Array of dual variables associated to the 1st set (rows of \p C)
 * \param[out] v Array of dual variables associated to the 2nd set (columns of \p C)
 * \param[in] init_type 0: no initialization, 1: classical (default)
 * \param[in] forb_assign If true, forbidden assignments are marked with negative values in the cost matrix
 * \details A solution to the LSAP is computed with the primal-dual version of the Hungarian algorithm, as detailed in:
 * \li <em>R. Burkard, M. Dell'Amico and S. Martello. Assignment Problems. SIAM, 2009</em>
 * \li <em>E.L. Lawler. Combinatorial Optimization: Networks and Matroids. Holt, Rinehart and Winston, New York, 1976</em>
 *
 * This version updates dual variables \c u and \c v instead of the cost matrix \c C, and at each iteration, the current matching is augmented by growing only one Hungarian tree until an augmenting path is found. Our implementation uses a Bread-First-like strategy to construct the tree, according to a FIFO strategy to select the next element at each iteration of the growing.
 *
 * Complexities:
 * \li O(min{n,m}²max{n,m}) in time (worst-case)
 * \li O(nm) in space
 *
 * \remark
 * Template \p DT allows to compute a solution with integer or floating-point values. Note that rounding errors may occur with floating point values when dual variables are updated but this does not affect the overall process. These errors can be observed when the reduced cost matrix is computed from the dual variables.
 */
// --------------------------------------------------------------------------------
template <class DT, typename IT>
void hungarianLSAP(const DT *C, const IT &n, const IT &m, IT *rho, DT *u, DT *v, IT *varrho = NULL, unsigned short init_type = 1, bool forb_assign = false)
{
  if (n == m) hungarianSquareLSAP<DT,IT>(C,n,rho,u,v,varrho,init_type,forb_assign);
  else hungarianRectLSAP<DT,IT>(C,n,m,rho,varrho,u,v,init_type,forb_assign);
}

// -----------------------------------------------------------
#endif
