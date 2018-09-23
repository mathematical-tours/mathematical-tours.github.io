function [bb] = cdtbal1(pp,ee)
%CDTBAL1 compute the circumballs associated with a 1-simplex
%triangulation embedded in R^2.
%   [BB] = TRIBAL1(PP,EE) returns the circumscribing balls
%   associated with the 1-simplexes in [PP,TT], such that BB 
%   = [XC,YC,RC.^2].

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 24/03/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ...
        ~isnumeric(ee) )
        error('cdtbal1:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(ee) ~= +2)
        error('cdtbal1:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(ee,2) < +2)
        error('cdtbal1:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    bb = zeros(size(ee,1),3);
    
    bb(:,1:2) = (pp(ee(:,1),:)+pp(ee(:,2),:))*.50 ;
    bb(:,  3) = ...
      sum((pp(ee(:,1),:)-pp(ee(:,2),:)).^2,2)*.25 ;
    
end



