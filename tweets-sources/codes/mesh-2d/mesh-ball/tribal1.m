function [bb] = tribal1(pp,ee)
%TRIBAL1 compute the circumballs associated with a 1-simplex
%triangulation embedded in R^2 or R^3.
%   [BB] = TRIBAL1(PP,EE) returns the circumscribing balls
%   associated with the edge segments in [PP,EE], such that 
%   BB = [XC,YC,RC.^2].

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 02/05/2018

    bb = pwrbal1(pp,zeros(size(pp,1),1),ee) ;

end



