% -----------------------------------------------------------
% file: compile_mex.m
% -----------------------------------------------------------
% authors: Sebastien Bougleux and Évariste Daller
% institution: Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC UMR 6072
% -----------------------------------------------------------
% Execute this file in matlab to compile matlab functions for LSAP and LSAPE
% -----------------------------------------------------------
% This file is part of LSAPE.
% LSAPE is free software: you can redistribute it and/or modify
% it under the terms of the CeCILL-C License. See README for more
% details.
% -----------------------------------------------------------

disp('Compiling tools for LSAP and LSAPE ...');
files =  { 'hungarianLSAP.cpp', ...
       	   'hungarianLSAPE.cpp', ...
           'greedyLSAP.cpp'
};

str = 'mex -O -I../include ';
flags = ' ';%' -std=c++11';

for i=1:length(files)
    eval([str files{i} flags ' ']);
end
