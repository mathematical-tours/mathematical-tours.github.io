function Ctrans = transformCost(C,u,v,neg)
%% 
% Ctrans = transformCost(C,u,v,neg)
%
% Compute Ctrans = C - (u1' + 1v)
% 
% -----------------------------------------------------------
% author: Sebastien Bougleux
% institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC
% ----------------------------------------------------------- 
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% ----------------------------------------------------------- 

    [n,m] = size(C);
    Ctrans = C - (repmat(u,1,m) + repmat(v,n,1));
    
    if exist('neg','var') && neg
        idxs = C < 0;
        Ctrans(idxs) = C(idxs);
    end
    
end
