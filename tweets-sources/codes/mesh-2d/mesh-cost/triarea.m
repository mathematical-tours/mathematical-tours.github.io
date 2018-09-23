function [area] = triarea(pp,tt)
%TRIAREA calc. triangle areas for a 2-simplex triangulation 
%embedded in the two-dimensional plane.
%   [AREA] = TRIAREA(VERT,TRIA) returns the signed triangle
%   areas, where AREA is a T-by-1 vector, VERT is a V-by-2 
%   array of XY coordinates, and TRIA is a T-by-3 array of
%   vertex indexing, where each row defines a triangle, such 
%   that VERT(TRIA(II,1),:), VERT(TRIA(II,2),:) and VERT(
%   TRIA(II,3),:) are the coordinates of the II-TH triangle.
%
%   See also TRISCR2, TRIANG2, TRIBAL2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 17/01/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(tt) )
        error('triarea:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(tt) ~= +2 )
        error('triarea:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(tt,2) < +3 )
        error('triarea:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nnod = size(pp,1) ;

%---------------------------------------------- basic checks
    if (min(min(tt(:,1:3))) < +1 || ...
            max(max(tt(:,1:3))) > nnod )
        error('triarea:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%--------------------------------------- compute signed area
    ev12 = pp(tt(:,2),:)-pp(tt(:,1),:) ;
    ev13 = pp(tt(:,3),:)-pp(tt(:,1),:) ;

    switch (size(pp,2))
        case +2
           
        area = ev12(:,1).*ev13(:,2) ...
             - ev12(:,2).*ev13(:,1) ;
        area = 0.5 * area;
            
        case +3
            
        avec = cross(ev12,ev13);
        area = sqrt(sum(avec.^2,2)) ;
        area = 0.5 * area;    
        
        otherwise
        error('Unsupported dimension.') ;
    end

end



