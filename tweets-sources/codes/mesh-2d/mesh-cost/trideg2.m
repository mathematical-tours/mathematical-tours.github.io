function [vdeg] = trideg2(pp,tt)
%TRIDEG2 calc. topological degree for vertices in a 2-simpl-
%ex triangulation.
%   [VDEG] = TRIDEG2(VERT,TRIA) returns the no. of triangles 
%   incident to each vertex. VDEG is a V-by-1 array of vert-
%   ex degrees, VERT is a V-by-D array of XY coordinates, 
%   and TRIA is a T-by-3 array of vertex indexing, where 
%   each row defines a triangle, such that 
%   VERT(TRIA(II,1),:), VERT(TRIA(II,2),:) and 
%   VERT(TRIA(II,3),:) are the coordinates of the II-TH tri-
%   angle.
%
%   See also TRISCR2, TRIVOL2, TRIANG2, TRIBAL2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 29/07/2018

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(tt) )
        error('trideg2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(tt) ~= +2 )
        error('trideg2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2) < +2 || size(tt,2) < +3 )
        error('trideg2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(pp,1) ;
    ntri = size(tt,1) ;

%---------------------------------------------- basic checks
    if (min(min(tt(:,1:3))) < +1 || ...
            max(max(tt(:,1:3))) > nvrt )
        error('trideg2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%------------------------------------- compute vertex degree
    vdeg = sum(sparse( ...
        tt(:,1:3),repmat( ...
            (1:ntri)',1,3),+1,nvrt,ntri),2) ;
            
    vdeg = full(vdeg);
    
end



