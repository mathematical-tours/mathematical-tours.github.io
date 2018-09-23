function [vert,tria,hlfs] = lfshfn2(varargin)
%LFSHFN2 calc. a discrete "local-feature-size" estimate for
%a polygonal domain embedded in R^2.
%   [VERT,TRIA,HFUN] = LFSHFN2(NODE,EDGE) returns the trian-
%   gulated "feature-size" estimate for the polygonal region 
%   {NODE,EDGE}. NODE is an N-by-2 array of polygonal verti-
%   ces and EDGE is an E-by-2 array of edge indexing. Each
%   row in EDGE represents an edge of the polygon, such that
%   NODE(EDGE(JJ,1),:) and NODE(EDGE(JJ,2),:) are the coord-
%   inates of the endpoints of the JJ-TH edge. If the argum-
%   ent EDGE is omitted it assumed that the vertices in NODE
%   are connected in ascending order.
%
%   [...] = LFSHFN2(NODE,EDGE,PART) computes a size-estimate
%   for a multiply-connected geometry. PART is a cell-array 
%   of polygonal "parts", where each element PART{KK} is an 
%   array of edge indices defining a given polygonal region. 
%   EDGE(PART{KK}, :) is the set of edges in the KK-TH part.
%
%   VERT is a V-by-2 array of XY coordinates, TRIA is a T-by
%   -3 array of triangles and HFUN is a V-by-1 array of mesh
%   -size values. Each row of TRIA defines a triangle, such 
%   that VERT(TRIA(II,1),:), VERT(TRIA(II,2),:) and VERT(
%   TRIA(II,3),:) are the coordinates of the II-TH triangle.
%
%   See also TRIHFN2, LIMHFN2, IDXTRI2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 24/01/2017

%---------------------------------------------- extract args
    node = []; PSLG = []; part = {}; opts = [] ; 
    
    if (nargin>=+1), node = varargin{1}; end
    if (nargin>=+2), PSLG = varargin{2}; end
    if (nargin>=+3), part = varargin{3}; end
    if (nargin>=+4), opts = varargin{4}; end
   
%------------------------------ build coarse background grid
   [opts] = makeopt(opts);
   
   [vert,conn,tria,tnum] = ...
        refine2(node,PSLG,part,opts) ;

%------------------------------ estimate local-feature-size!
    hlfs = ...
    +inf * ones(size(vert,1),1) ;

%------------------------------ calc. LFS based on edge-len.
    evec = vert(conn(:,2),:) - ...
           vert(conn(:,1),:) ;
    elen = sqrt(sum(evec.^2,2)) ;
    hlen = elen * +1. ;
   
    for epos = +1 : size(conn,+1)
    
        ivrt = conn(epos,1) ;
        jvrt = conn(epos,2) ;  
        
        hlfs(ivrt) = min( ...
          hlfs(ivrt),hlen(epos)) ;
        hlfs(jvrt) = min( ...
          hlfs(jvrt),hlen(epos)) ;
    
    end

%------------------------------ push gradient limits on HFUN   
    DHDX = opts.dhdx;
   
    hlfs = ...
        limhfn2(vert,tria,hlfs,DHDX) ;

end

function [opts] = makeopt(opts)
%MAKEOPT setup the options structure for LFSHFN2.

    if (~isfield(opts,'kind'))
        opts.kind = 'delaunay';
    else
    if (~strcmpi(opts.kind,'delfront') && ...
        ~strcmpi(opts.kind,'delaunay') )
        error( ...
        'lfshfn2:invalidOption','Invalid refinement KIND.'); 
    end
    end
    
    if (~isfield(opts,'rho2'))
        opts.rho2 = sqrt(+2.) ;
    else
    if (~isnumeric(opts.rho2))
        error('lfshfn2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.rho2)~= +1)
        error('lfshfn2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.rho2 < +1.)
        error('lfshfn2:invalidOptionValues', ...
            'Invalid OPT.RHO2 selection.') ;
    end
    end
    
    if (~isfield(opts,'dhdx'))
        opts.dhdx = +0.2500 ;
    else
    if (~isnumeric(opts.dhdx))
        error('lfshfn2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.dhdx)~= +1)
        error('lfshfn2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.dhdx <= 0.)
        error('lfshfn2:invalidOptionValues', ...
            'Invalid OPT.DHDX selection.') ;
    end
    end

end



