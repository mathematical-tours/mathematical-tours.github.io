% hungarianLSAP.m Help file for hungarianLSAP MEX-file.
%  hungarianLSAP.c - Compute a solution to (symmetric or asymmetric) LSAP with Hungarian algorithm
% 
%             [rho] = hungarianLSAP(C,init_type,forb)
%      [rho,varrho] = hungarianLSAP(C,init_type,forb)
%         [rho,u,v] = hungarianLSAP(C,init_type,forb)
%  [rho,varrho,u,v] = hungarianLSAP(C,init_type,forb)
%
%  Given a nxm cost matrix C (integer ou floatting values) it computes:
%
%  - a solution rho to the LSAP, i.e. a one-to-one mapping from the rows of C to its columns,
%    and the one-to-to mapping from the columns to the rows
%
%    rho is represented by a nx1 matrix so that rho(i)=j means that i is assigned to j
%    varrho is a 1xm matrix so that varrho(j)=i means that j is assigned to i
%    for the asymmetric LSAP (n and m are not equal):
%      if n<m:
%         rho is an injection from {1,...,n} to {1,...,m}
%         m-n columns are not assigned, which is represented by varrho(j)=-1
%      if n>m
%         varrho is an injection from {1,...,m} to {1,...,n}
%         n-m rows are not assigned, which is represented by rho(i)=-1
%
%  - a solution (u,v) to its dual problem (labeling problem)
%
%  optional init_type:
%    0: no initialization (u=v=0)
%    1 (default): classical initialization (u=min(C) on rows, v=min(C-u) on columns)
%
%  optional boolean parameter forb:
%    true  -> forbidden assignments are represented by negative cost values
%    false -> no forbidden assignments (by default) 
%
%  This is a MEX-file for MATLAB.
%  This file is part of LSAPE.
%  LSAPE is free software: you can redistribute it and/or modify it
%  under the terms of the CeCILL-C License. See README for more details.
%
%     Copyright 2015-2017
%      authors: Sebastien Bougleux
%  institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, France
%   last modif: July 5 2017
%
