/* This file is a matlab mex function for computing the transportation cost 
 * between two vectors given a cost matrix.
 *
 * It was written by Antoine Rolet (2014) and mainly consists of a wrapper
 * of the code written by Nicolas Bonneel available on this page
 *          http://people.seas.harvard.edu/~nbonneel/FastTransport/
 *
 * Please give relevant credit to the original author (Nicolas Bonneel) if
 * you use this code for a publication.
 *
 */
#include "mex.h"
#include <iostream>
#include <vector>

#include "network_simplex_simple.h"


using namespace lemon;
typedef unsigned int node_id_type;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int n, m, i, n1, n2,cur;
  double *X, *Y, *D, *demd, *optimal, max,max_iter;
  const mwSize *dimx, *dimy, *dimd, *dimw;
	
    
  typedef FullBipartiteDigraph Digraph;
  DIGRAPH_TYPEDEFS(FullBipartiteDigraph);    
    
  // Check if the number of input and output argument is correct
    
  if (nrhs < 3)
    mexErrMsgTxt("There should be at least 3 input arguments");
  if (nlhs > 2)
    mexErrMsgTxt("Too many output argument");
    
  // Get the input arguments
    
  dimx = mxGetDimensions(prhs[0]);
  dimy = mxGetDimensions(prhs[1]);
  dimd = mxGetDimensions(prhs[2]);
  X = mxGetPr(prhs[0]);
  Y = mxGetPr(prhs[1]);
  D = mxGetPr(prhs[2]);

  if(nrhs==4)
    max_iter=*mxGetPr(prhs[3]);
  else
    max_iter=-1;
    
  // Check dimensions
    
  if (dimx[1] != 1 || dimy[1] != 1)
    mexErrMsgTxt("First two arguments should be column vectors");
  n1 = dimd[0];
  n2 = dimd[1];
  if (dimx[0] != n1 || dimy[0] != n2)
    mexErrMsgTxt(
		 "Ground metric matrix dimensions should be the length of the two first arguments");
  
  // Get the number of non zero coordinates for r and c
    n=0;
    for (node_id_type i=0; i<n1; i++) {
        double val=*(X+i);
        if (val>0) {
            n++;
        }
    }
    m=0;
    for (node_id_type i=0; i<n2; i++) {
        double val=*(Y+i);
        if (val>0) {
            m++;
        }
    }
    
    
    // Define the graph
    
    std::vector<int> indI(n), indJ(m);
    std::vector<double> weights1(n), weights2(m);
    Digraph di(n, m);
    NetworkSimplexSimple<Digraph,double,double, node_id_type> net(di, true, n+m, n*m,max_iter);
    
    // Set supply and demand, don't account for 0 values (faster)
    
    max=0;
    cur=0;
    for (node_id_type i=0; i<n1; i++) {
        double val=*(X+i);
        if (val>0) {
            weights1[ di.nodeFromId(cur) ] = val;
            max+=val;
            indI[cur++]=i;
        }
    }
    
    // Demand is actually negative supply...
    
    max=0;
    cur=0;
    for (node_id_type i=0; i<n2; i++) {
        double val=*(Y+i);
        if (val>0) {
            weights2[ di.nodeFromId(cur) ] = -val;
            indJ[cur++]=i;
            
            max-=val;
        }
    }
    
    
    net.supplyMap(&weights1[0], n, &weights2[0], m);
    
    // Set the cost of each edge
    max=0;
    for (node_id_type i=0; i<n; i++) {
        for (node_id_type j=0; j<m; j++) {
            double val=*(D+indI[i]+indJ[j]*n1);
            net.setCost(di.arcFromId(i*m+j), val);
            if (val>max) {
                max=val;
            }
        }
    }
    
    
    // Solve the problem with the network simplex algorithm
    
    int ret=net.run();
    if (ret!=(int)net.OPTIMAL) {
        if (ret==(int)net.INFEASIBLE) {
            mexErrMsgTxt("Infeasible problem");
        }
        if (ret==(int)net.UNBOUNDED) {
            mexErrMsgTxt("Unbounded problem");
        }
    }
    
    // Define the output arguments
    
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	demd = mxGetPr(plhs[0]);
    
	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(n1, n2, mxREAL);
		optimal = mxGetPr(plhs[1]);
	}
    
    // Set the output values
	
	if (nlhs > 1) {
		for (node_id_type i=0; i<n; i++) {
			for (node_id_type j=0; j<m; j++)
			{
				*(optimal+indI[i]+indJ[j]*n1) = net.flow(di.arcFromId(i*m+j));
			}
		}
	}
    
    *demd=net.totalCost();
    
}
