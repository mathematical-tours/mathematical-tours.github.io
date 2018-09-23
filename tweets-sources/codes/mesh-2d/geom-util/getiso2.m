function [node,edge] = ...
        getiso2(xpos,ypos,zdat,ilev,filt)
%GETISO2 extract an iso-contour from a structured two-dimen-
%sional data-set. 
%   [NODE,EDGE] = GETISO2(XPOS,YPOS,ZFUN,ZLEV) returns the 
%   contour ZFUN(XPOS,YPOS) = ZLEV as a PSLG, by post-proce-
%   ssing the output of the CONTOUR function. The arguments
%   XPOS, YPOS and ZFUN must all be N-by-M arrays, and ZLEV 
%   a scalar contouring value.
%   
%   See also GETNAN2, FIXGEO2, BFSGEO2, REFINE2

%-----------------------------------------------------------
%   Darren Engwirda : 2018 --
%   Email           : de2363@columbia.edu
%   Last updated    : 25/03/2018
%-----------------------------------------------------------

    if (nargin < +5), filt = +0. ; end
    
%---------------------------------------------- basic checks    
    if ( ~isnumeric(xpos) || ...
         ~isnumeric(ypos) || ...
         ~isnumeric(zdat) || ...
         ~isnumeric(ilev) || ...
         ~isnumeric(filt) )
        error('getiso2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(xpos) ~= +2 || ...
        ndims(ypos) ~= +2 || ...
        ndims(zdat) ~= +2 )
        error('getiso2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (isvector(xpos))
        xnum = length(xpos);
    else
        xnum = size(xpos,2);
    end
    
    if (isvector(ypos))
        ynum = length(ypos);
    else
        ynum = size(ypos,1);
    end
    
    if (xnum ~= size(zdat,2) || ...
        ynum ~= size(zdat,1) || ...
        numel(ilev) ~= +1 || ...
        numel(filt) ~= +1 )
        error('getiso2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

%------------------------------------ compute the isocontour
    cmat = contourc( ...
        xpos,ypos,zdat,[ilev,ilev]) ;

%------------------------------------ "walk" contour segment
    node = [] ; edge = [] ; ipos = +1 ;

    while (ipos < size(cmat,2))
       
        numc = cmat(2,ipos);
        ppts =[cmat(1,ipos+1:ipos+numc)', ...
               cmat(2,ipos+1:ipos+numc)'
              ] ;

        pmin = min(ppts,[],1);
        pmax = max(ppts,[],1);
        
        pdel = pmax - pmin ;
        
        if (min(pdel)>=filt)
        
            if all(ppts(1,:) == ppts(end,:))

    %-------------------------------- closed - back to start
            enew = ...
           [(1:numc-1)',(2:numc-0)'; numc,1] ;

            else
            
    %-------------------------------- open - dangling endpts
            enew = ...
           [(1:numc-1)',(2:numc-0)'] ;

            end

            enew = ...
            enew + size(node,1);

            node = [node; ppts];
            edge = [edge; enew];

        end

        ipos = ipos + numc + 1 ;
        
    end

end



