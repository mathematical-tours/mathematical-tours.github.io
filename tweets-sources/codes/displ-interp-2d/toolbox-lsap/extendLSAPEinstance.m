function [Cext,val] = extendLSAPEinstance(C,val)
% Extend a cost matrix for error-correcting matching to its square version
%
%  [Cext,val] = extendCostMtx(C,val)
%
%  inputs:
%  C (n+1)x(m+1) cost matrix of an error-correcting bipartite graph
%  val the value used to fill 2nd and 3rd block of Cext (optional)
%  val = max(max(C)) + 1 if not given
%
%  outputs:
%  Cext (n+m)x(m+n) extended cost matrix
%  valcpt (optional) is the value val used to fill 2nd and 3rd block of
%  Cext
%  
%         | c(1,1)   ...  c(1,m-1)  | c(1,m) val   ...      val   |
%         |   .              .      |  val  c(2,m) val ...   .    | 
%         |   .              .      |   .                   val   |
%         | c(n-1,1) ... c(n-1,m-1) |  val     ...   val c(n-1,m) |
%  Cext =  -------------------------------------------------------
%         | c(n,1) val  ...    val  |                             |
%         |  val c(n,2) val ..  .   |           0_{m,n}           |
%         |   .                val  |                             |
%         |  val  ...  val c(n,m-1) |                             |
%
%
%   author: Sebastien Bougleux
%   institution: Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC
% -----------------------------------------------------------
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% -----------------------------------------------------------
    
    if ~exist('val','var') || nargin == 1
        val = max(max(C)) + 1;
    end
    
    D = val .* cast(ones(size(C,1)-1),class(C));
    E = val .* cast(ones(size(C,2)-1),class(C));
    D(logical(eye(size(D)))) = C(1:end-1,end);
    E(logical(eye(size(E)))) = C(end,1:end-1);
    Cext = [[C(1:end-1,1:end-1),D];[E,zeros(size(C,2)-1,size(C,1)-1)]];

end
