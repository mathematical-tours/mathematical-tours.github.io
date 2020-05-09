% Legendre-Fenchel Computation Package
% Version 1.0  19-Sept.-97.
% Copyright (c) 1997 by Yves Lucet.
%
%
% Legendre-Fenchel computation routines.
%
%   LFt       - Compute the Legendre Transform of a function with the LLT
%   LLTd      - Compute the Legendre Transform of discrete data
%   LLTd2D    - Compute the Legendre Transform of discrete data
%		corresponding to a bivariate function
%   LFtdirect - Compute the Legendre Transform of a function directly
%   LLTdirect - Compute the Legendre Transform of discrete data directly
%
%
% Convex computation routines.
%
%   bb        - Compute the convex hull of discrete data
%   infconvd  - Compute the inf-convolution of convex discrete data
%
%
% Auxiliary routines required by some of the above routines.
%
%   fusion    - Switch to fusionma. It can easily be modified to
%		switch to fusionca.
%   fusionma  - Use matlab syntax to quickly merge two increasing sequences.
%   fusionca  - Classical programming to merge two increasing sequences.
%
