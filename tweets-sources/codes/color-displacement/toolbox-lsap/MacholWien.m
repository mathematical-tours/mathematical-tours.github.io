function C = MacholWien( n,m,c )
%
%   C = MacholWien(n,m,c)
%
%   compute the nxm matrix C sot that c(i,j)=c*i*j
%   aka Machol-Wien instance for the LSAP
%   known to be a worst-case time complexity for several algorithms
%
% authors: Sebastien Bougleux
% institution: Normandie Univ, CNRS - UNICAEN - ENSICAEN, GREYC, France
% ----------------------------------------------------------- 
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% ----------------------------------------------------------- 

    if nargin < 3
        c = 1;
    end
        
    C = repmat(1:m,n,1);
    C = bsxfun(@times,C,c.*(1:n)');

end

