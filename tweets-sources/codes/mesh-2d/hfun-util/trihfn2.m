function [hval] = trihfn2(test,vert,tria,tree,hfun)
%TRIHFN2 evaluate a discrete mesh-size function defined on a 
%2-simplex triangulation embedded in R^2.
%   [HVAL] = TRIHFN2(TEST,VERT,TRIA,TREE,HFUN) returns an
%   interpolation of the discrete mesh-size function HFUN to
%   the queries points TEST. HVAL is a Q-by-1 array of func-
%   tion values associated with the Q-by-2 array of XY coor-
%   dinates TEST. {VERT,TRIA,HFUN} is a discrete mesh-size
%   function, where VERT is a V-by-2 array of XY coordinates
%   TRIA is a T-by-3 array of triangles and HFUN is a V-by-1
%   array of mesh-size values. Each row of TRIA defines a 
%   triangle, such that VERT(TRIA(II,1),:), 
%   VERT(TRIA(II,2),:) and VERT(TRIA(II,3),:) are the coord-
%   inates of the II-TH triangle. TREE is a spatial indexing
%   structure for {VERT,TRIA}, as returned by IDXTRI2.
%
%   See also LFSHFN2, LIMHFN2, IDXTRI2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 07/04/2017

%---------------------------------------------- basic checks    
    if ( ~isnumeric(test) || ...
         ~isnumeric(vert) || ...
         ~isnumeric(tria) || ...
         ~isnumeric(hfun) || ...
         ~isstruct (tree) )
        error('trihfn2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(test) ~= +2 || ...
        ndims(vert) ~= +2 || ...
        ndims(tria) ~= +2 || ...
        ndims(hfun) ~= +2 )
        error('trihfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(test,2)~= +2 || ...
        size(vert,2)~= +2 || ...
        size(tria,2) < +3 || ...
        size(hfun,2)~= +1 || ...
        size(vert,1)~= size(hfun,1) )
        error('trihfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(vert,1) ;

%---------------------------------------------- basic checks
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nvrt )
        error('trihfn2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%-------------------------------------- test-to-tria queries
   [tp,tj] = ...
    findtria (vert,tria,test,tree);
   
    if (isempty(tp))
        in = false(size(test,1),1);
        ti = [];
    else
        in = tp(:,1) > +0 ;
        ti = tj(tp(in,+1));
    end

%-------------------------------------- calc. linear interp.    
    hval = max(hfun) * ones(size(test,1),1) ;

    if (any(in))

    d1 = test(in,:) - vert(tria(ti,1),:);
    d2 = test(in,:) - vert(tria(ti,2),:);
    d3 = test(in,:) - vert(tria(ti,3),:);
    
    a3 = abs(d1(:,1) .* d2(:,2) ...
           - d1(:,2) .* d2(:,1) ) ;
    a2 = abs(d1(:,1) .* d3(:,2) ...
           - d1(:,2) .* d3(:,1) ) ;
    a1 = abs(d3(:,1) .* d2(:,2) ...
           - d3(:,2) .* d2(:,1) ) ;
    
    hval(in) = a1.*hfun(tria(ti,1)) ...
             + a2.*hfun(tria(ti,2)) ...
             + a3.*hfun(tria(ti,3)) ;
    
    hval(in) = hval(in)./(a1+a2+a3) ;

    end

end



