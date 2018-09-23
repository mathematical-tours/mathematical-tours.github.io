function [hrel] = relhfn2(vert,tria,hvrt)
%RELHFN2 calc. relative edge-length for a 2-simplex triangu-
%lation embedded in euclidean space.
%   [HREL] = RELHFN2(VERT,TRIA,HVRT) returns the relative 
%   edge-length, indicating conformance to the imposed mesh-
%   spacing constraints, where HREL is a E-by-1 vector, VERT 
%   is a V-by-D array of XY coordinates, TRIA is a T-by-3 
%   array of vertex indexing, and HVRT is a V-by-1 array of
%   mesh-spacing values associated with the mesh vertices.  
%
%   See also TRISCR2, TRIVOL2, TRIANG2, TRIBAL2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 10/08/2018

%---------------------------------------------- basic checks    
    if (~isnumeric(vert) || ~isnumeric(tria) ||   ...
        ~isnumeric(hvrt) )
        error('relhfn2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ndims(tria) ~= +2)
        error('relhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || size(tria,2) < +3)
        error('relhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(hvrt,2)~= +1 || size(hvrt,1) ~= size(vert,1) )
        error('relhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nnod = size(vert,1) ;

%---------------------------------------------- basic checks
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nnod )
        error('relhfn2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%----------------------------------- compute rel. mesh-sizes
    eset = unique2([ ...
        tria(:,[1,2]); tria(:,[2,3])
        tria(:,[3,1])
                   ]);
   
    evec = vert(eset(:,2),:)-vert(eset(:,1),:) ;
         
    elen = sqrt(sum(evec.^2,+2));

    hmid = hvrt(eset(:,2),:)+hvrt(eset(:,1),:) ;
    hmid = hmid * +0.50 ;
    hrel = elen ./ hmid ;

end


