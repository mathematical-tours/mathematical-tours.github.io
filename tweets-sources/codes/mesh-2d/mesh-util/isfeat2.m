function [is,bv] = isfeat2(pp,ee,tt)
%ISFEAT2 return "feature" status for the triangles in a two-
%dimensional constrained triangulation.
%   [STAT] = ISFEAT2(VERT,EDGE,TRIA) returns STAT(II) = TRUE
%   for any triangle including a sufficiently "sharp" angle 
%   located at the apex of any two constrained edges. Sharp
%   features have angles of greater-than ACOS(+0.8) degrees.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 27/01/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(ee) || ...
        ~isnumeric(tt) )
        error('isfeat2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(ee) ~= +2 || ...
        ndims(tt) ~= +2 )
        error('isfeat2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(ee,2) < +5 || ...
        size(tt,2) < +6 )
        error('isfeat2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nnod = size(pp,1) ; nedg = size(ee,1) ;
    ntri = size(tt,1) ;

%---------------------------------------------- basic checks
    if (min(min(tt(:,1:3))) < +1 || ...
        max(max(tt(:,1:3))) > nnod)
        error('isfeat2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end
    if (min(min(tt(:,4:6))) < +1 || ...
        max(max(tt(:,4:6))) > nedg)
        error('isfeat2:invalidInputs', ...
            'Invalid TRIA input array.') ;   
    end
    
    if (min(min(ee(:,1:2))) < +1 || ...
        max(max(ee(:,1:2))) > nnod)
        error('isfeat2:invalidInputs', ...
            'Invalid EDGE input array.') ;   
    end
    if (min(min(ee(:,3:4))) < +0 || ...
        max(max(ee(:,3:4))) > ntri)
        error('isfeat2:invalidInputs', ...
            'Invalid EDGE input array.') ;   
    end

%----------------------------------------- compute "feature"
    is = false(size(tt,1),1);
    bv = false(size(tt,1),3);

    EI = [3, 1, 2] ;
    EJ = [1, 2, 3] ;    
    NI = [3, 1, 2] ;
    NJ = [1, 2, 3] ;
    NK = [2, 3, 1] ;

    for ii = +1 : +3
    
    %------------------------------------- common edge index
        ei = tt( :,EI(ii)+3);
        ej = tt( :,EJ(ii)+3);
    
    %------------------------------------- is boundary edge? 
        bi = ee(ei,5) >= +1 ;
        bj = ee(ej,5) >= +1 ;
        
        ok = bi & bj ;
        
        if (~any(ok)), continue; end
        
        ni = tt(ok,NI(ii)+0);
        nj = tt(ok,NJ(ii)+0);
        nk = tt(ok,NK(ii)+0);
        
    %------------------------------------- adj. edge vectors
        vi = pp(ni,:)-pp(nj,:) ;
        vj = pp(nk,:)-pp(nj,:) ;

    %------------------------------------- adj. edge lengths
        li = sqrt(sum(vi.^2,2));
        lj = sqrt(sum(vj.^2,2));

        ll = li .* lj ;
            
    %------------------------------------- adj. dot-product!
        aa = sum(vi.*vj,2)./ll ;
    
        bv(ok,ii) = aa >= +.80 ;
    
        is(ok) = ...
            is(ok) | bv(ok,ii) ;
 
    end

end



