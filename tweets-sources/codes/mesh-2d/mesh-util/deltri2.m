function [vert,conn,tria,tnum] = deltri2(varargin)
%DELTRI2 compute a constrained 2-simplex Delaunay triangula-
%tion in the two-dimensional plane.
%   [VERT,CONN,TRIA,TNUM]=DELTRI2(VERT,CONN,NODE,PSLG,PART)
%   computes the Delaunay trianguation {VERT,TRIA}, the con-
%   straints CONN, and the "inside" status vector TNUM. VERT
%   is an V-by-2 array of XY coordinates to be triangulated,
%   TRIA is a T-by-3 array of vertex indexing, where each 
%   row defines a triangle, such that VERT(TRIA(II,1),:), 
%   VERT(TRIA(II,2),:) and VERT(TRIA(II,3),:) are the coord-
%   inates of the II-TH triangle. CONN is a C-by-2 array of
%   constraining edges, where each row defines an edge, as
%   per TRIA. The additional arguments NODE,PSLG and PART 
%   define a (mutliply-connected) polygonal region, where 
%   NODE is an N-by-2 array of vertices and PSLG is a P-by-2 
%   array of edges (a piecewise-straight-line-graph), where 
%   each row defines an edge as a pair of indices into NODE.
%   PART is a cell-array of polygonal "parts", where each
%   element PART{KK} is an array of edge indices defining a
%   polygonal region. PSLG(PART{KK},:) is the set of edges
%   in the KK-TH part. TNUM is a T-by-1 array of part index-
%   ing, such that TNUM(II) is the index of the part in whi-
%   ch the II-TH triangle resides.
%
%   See also DELAUNAYTRIANGULATION, DELAUNAYTRI, DELAUNAYN

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 08/07/2018

    vert = []; conn = []; node = []; PSLG = [];
    part = {}; kind = 'constrained'; 
    
%---------------------------------------------- extract args
    if (nargin>=+1), vert = varargin{1}; end
    if (nargin>=+2), conn = varargin{2}; end
    if (nargin>=+3), node = varargin{3}; end
    if (nargin>=+4), PSLG = varargin{4}; end
    if (nargin>=+5), part = varargin{5}; end
    if (nargin>=+6), kind = varargin{6}; end
    
%---------------------------------------------- basic checks    
    if (~isnumeric(vert) || ~isnumeric(conn) || ...
        ~isnumeric(node) || ~isnumeric(PSLG) || ...
        ~iscell   (part) || ~ischar   (kind) )
        error('deltri2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

    nvrt = size(vert,+1) ; nnod = size(node,+1) ;
    nedg = size(PSLG,+1) ;

%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ndims(conn) ~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || size(conn,2)~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (min([conn(:)])<+1 || max([conn(:)])>nvrt)
        error('deltri2:invalidInputs', ...
            'Invalid CONN input array.') ;
    end
    
%---------------------------------------------- basic checks
    if (nargin >= +3)

    if (ndims(node) ~= +2 || ndims(PSLG) ~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(node,2)~= +2 || size(PSLG,2)~= +2)
        error('deltri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    if (min([PSLG(:)])<+1 || max([PSLG(:)])>nnod)
        error('deltri2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    pmin = cellfun(@min,part);
    pmax = cellfun(@max,part);
    
    if (min([pmin(:)])<+1 || max([pmax(:)])>nedg)
        error('deltri2:invalidInputs', ...
            'Invalid PART input array.') ;
    end

    end

%------------------------------------ compute Delaunay tria.
    switch (lower(kind))
    case 'constrained'

        if (exist( ...
        'delaunayTriangulation') == +2 )
    %-------------------------------- use class if available    
        dtri = ...
        delaunayTriangulation(vert,conn) ;
        vert = dtri.Points;
        conn = dtri.Constraints;
        tria = dtri.ConnectivityList;
        else
        if (exist('DelaunayTri') == +2 )
    %-------------------------------- use class if available
        dtri = DelaunayTri   (vert,conn) ;
        vert = dtri.X;
        conn = dtri.Constraints;
        tria = dtri.Triangulation;
        else
    %-------------------------------- *fall-back* onto qhull
       [vert,conn,tria] ...
                = cfmtri2(vert,conn) ;
        end
        end
    
    case 'conforming'   
        
    %-------------------------------- "conforming" delaunay!
       [vert,conn,tria] ...
                = cfmtri2(vert,conn) ;
        
    otherwise
        error('deltri2:invalidInputs', ...
            'Invalid KIND selection.') ;
    
    end
   
%------------------------------------ calc. "inside" status!    
    tnum = zeros(size(tria,+1),+1) ;
    
    if (nargin >= +3)
    
    tmid = vert(tria(:,1),:) ...
         + vert(tria(:,2),:) ...
         + vert(tria(:,3),:) ;
    tmid = tmid / +3.0;

    for ppos = 1 : length(part)

       [stat] = inpoly2( ...
            tmid,node  , ...
            PSLG(part{ppos},:))  ;
 
        tnum(stat)  = ppos ;
        
    end
    
%------------------------------------ keep "interior" tria's  
    tria = tria(tnum>+0,:) ;
    tnum = tnum(tnum>+0,:) ;
    
    end
    
%------------------------------------ flip for correct signs    
    area = triarea(vert,tria) ;

    tria(area<0.,:) = ...
        tria(area<0.,[1,3,2]) ;

end



