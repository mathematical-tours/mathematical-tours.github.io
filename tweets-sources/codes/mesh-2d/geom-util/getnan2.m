function [node,edge] = getnan2(varargin)
%GETNAN2 parse a NaN delimited polygon into a PSLG.
%   [NODE,EDGE] = GETNAN2(NANS,FILT) converts a set of NaN 
%   delimited polygons to a PSLG representation. NANS is an 
%   D-by-2 array of coordinates, with polygon vertices spec-
%   ified in connsecutive order, and delimited by NaN values
%   where breaks between polygons occur. FILT is a length 2
%   vector of "small-feature" filter values. Small polygons
%   wiith axis-aligned extents less than FILT are stripped
%   from the output. NODE is an N-by-2 array of coordinates
%   and EDGE is an E-by-2 array of edge indexing. Each row 
%   in EDGE represents an edge of the polygon, such that
%   NODE(EDGE(JJ,1),:) and NODE(EDGE(JJ,2),:) are the coord-
%   inates of the endpoints of the JJ-TH edge. 
%   
%   See also FIXGEO2, BFSGEO2, REFINE2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 06/10/2017
%-----------------------------------------------------------

    data = [] ; filt = +0. ;

    if (nargin>=+1), data = varargin{1}; end
    if (nargin>=+2), filt = varargin{2}; end
    
%---------------------------------------------- basic checks    
    if ( ~isnumeric(data) || ...
         ~isnumeric(filt) )
        error('getnan2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(data) ~= +2 )
        error('getnan2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (size(data,2)~= +2 || ...
        size(filt,2)>= +2 )
        error('getnan2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end  
    
%---------------------------------- parse NaN delimited data
    nvec = find(isnan(data(:,1))) ;
    
    if (isempty(nvec))                    % no NaN's at all!
    nvec = [nvec ; size(data,1) ] ;
    end
    
    if (nvec(end)~=size(data,1) )         % append last poly
    nvec = [nvec ; size(data,1) ] ;
    end
  
    node = zeros(size(data,1), 2) ;
    edge = zeros(size(data,1), 2) ;
    
    next = +1; nout = +1; eout = +1 ;
    
    for npos = +1 : length(nvec)
        
        stop = nvec(npos) ;
        
        pnew = data(next:stop-1,1:2);
        
        pmin = min(pnew,[],1) ;
        pmax = max(pnew,[],1) ;
        
        pdel = pmax-pmin;
        
        if (~isempty(pnew))
        if (any(pdel>filt))            
%---------------------------------- push polygon onto output            
        nnew = size(pnew,1);
        
        enew = [(1:nnew-1)', ...
                (2:nnew-0)'; ...
                nnew, +1 ] ;
        enew = enew+nout-1 ;
        
        mnew = size(enew,1);
        
        pvec = nout:nout+nnew-1;
        evec = eout:eout+mnew-1;
        
        node(pvec,:) = pnew;
        edge(evec,:) = enew;

        nout = nout + nnew ;
        eout = eout + mnew ;
        
        end
        end
        
        next = stop + 1 ;
        
    end

    node = node(+1:nout-1,:) ;
    edge = edge(+1:eout-1,:) ;
    
end


