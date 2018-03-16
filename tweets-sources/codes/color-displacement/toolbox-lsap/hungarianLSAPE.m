% hungarianLSAPE.m Help file for hungarianLSAP MEX-file.
%  hungarianLSAPE.cpp - Compute a solution to (symmetric or asymmetric) LSAP with Hungarian algorithm
% 
%      [rho,varrho] = hungarianLSAPE(C,init_type,forb)
%         [rho,u,v] = hungarianLSAPE(C,init_type,forb)
%  [rho,varrho,u,v] = hungarianLSAPE(C,init_type,forb)
%
%  Given a nxm cost matrix C (integer ou floatting values)
%  with last column encoding removal costs and last row insertion costs, it computes:
%
%  - a solution rho to the LSAPE, i.e. a mapping from the rows of C to its columns,
%    and the mapping varrho from the columns to the rows
%
%    rho is represented by a (n-1)x1 matrix so that rho(i)=j means that i is assigned to j
%    varrho is a 1x(m-1) matrix so that varrho(j)=i means that j is assigned to i
%    rho(i)=m or varrho(j)=n encode assignments to the null element (removal or insertion)
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
