% -----------------------------------------------------------
% file: nbECGM.m
% -----------------------------------------------------------
% authors: Sebastien Bougleux (UNICAEN) and Luc Brun (ENSICAEN)
% institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072 
% ----------------------------------------------------------- 
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README file 
% for more details.
% ----------------------------------------------------------- 

function nb = nbECGM(nbU,nbV)

    nb = 0;
    
    for p=0:min(nbU,nbV)
        nb = nb + factorial(p) * nchoosek(nbU,p) * nchoosek(nbV,p);
    end

end
