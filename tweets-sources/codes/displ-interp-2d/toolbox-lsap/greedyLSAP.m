% greedyLSAP.m Help file for greedyLSAP MEX-file.
%  greedyLSAP.cpp - Compute an approximate solution to (symmetric or
%  asymmetric) LSAP with greedy algorithms, i.e. an assignment having a
%  lost cost
% 
%           [rho,cost] = greedyLSAP(C,greedy_type)
%    [rho,varrho,cost] = greedyLSAP(C,greedy_type)
%
%  Given a nxm cost matrix C (integer ou floatting values), with n>=m, 
%  it computes an assignment with low cost
%
%  rho is the assignment from the rows to the columns (nx1 matrix)
%  varrho is the assignment from the columns to the rows (1xm matrix)
%  cost is the cost of the assignment
%  
%  optional integer greedy_type:
%  0: Basic
%  1: Refined
%  2: Loss (default)
%  3: Basic sort (C++ std::sort)
%  4: Counting sort (integers only)
% 
%   supporting cost values: int16, int32, int64, single, double (default)
%
% 
%   This is a MEX-file for MATLAB.
%   
%   This file is part of LSAPE.
%   LSAPE is free software: you can redistribute it and/or modify it
%   under the terms of the CeCILL-C License. See README for more details.
% 
%      Copyright 2015-2017
%       authors: Sebastien Bougleux
%   institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, France
%    last modif: July 6 2017
%  
%   execute matlab file 'compile_mex.m' to compile this function and use it in matlab
