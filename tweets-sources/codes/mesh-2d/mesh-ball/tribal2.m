function [bb] = tribal2(pp,tt)
%TRIBAL2 compute the circumballs associated with a 2-simplex
%triangulation embedded in R^2 or R^3.
%   [BB] = TRIBAL2(PP,TT) returns the circumscribing balls
%   associated with the triangles in [PP,TT], such that BB = 
%   [XC,YC,RC.^2].

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 02/05/2018

    bb = pwrbal2(pp,zeros(size(pp,1),1),tt) ;

end



