function [ll,ei] = minlen2(pp,tt)
%MINLEN2 return the minimum length edge for each triangle in 
%a two-dimensional triangulation.
%   [ELEN,IMIN] = MINLEN2(VERT,TRIA) returns the minimum le-
%   ngth ELEN and local edge index IMIN for all triangles in
%   the triangulation {VERT,TRIA}.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 16/01/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(tt) )
        error('minlen2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(tt) ~= +2 )
        error('minlen2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(tt,2) < +3 )
        error('minlen2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nnod = size(pp,1) ;

%---------------------------------------------- basic checks
    if (min(min(tt(:,1:3))) < +1 || ...
            max(max(tt(:,1:3))) > nnod )
        error('minlen2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end
    
%------------------------------------------ compute edge-len
    l1 = sum((pp(tt(:,2),:) ...
             -pp(tt(:,1),:)).^2,2);
    l2 = sum((pp(tt(:,3),:) ...
             -pp(tt(:,2),:)).^2,2);
    l3 = sum((pp(tt(:,1),:) ...
             -pp(tt(:,3),:)).^2,2);

%------------------------------------------ compute min.-len
   [ll,ei] = min([l1,l2,l3],[],+2);

end



