function [hfun] = limhfn2(vert,tria,hfun,dhdx)
%LIMHFN2 impose gradient limits on a discrete mesh-size fun-
%ction defined over a 2-simplex triangulation.
%   [HFUN] = LIMHFN2(VERT,TRIA,HFUN,DHDX) returns a "gradie-
%   nt-limited" function HFUN, defined over a triangulation 
%   {VERT,TRIA}. HFUN is a T-by-1 vector of function values, 
%   VERT is a V-by-2 array of XY coordinates and TRIA is a 
%   T-by-3 array of triangles. Each row of TRIA 
%   defines a triangle, such that VERT(TRIA(II,1),:), VERT(
%   TRIA(II,2),:) and VERT(TRIA(II,3),:) are the coordinates 
%   of the II-TH triangle. DHDX is a scalar gradient-limit.
%   HFUN is "limited" to control variation over the elements
%   in the triangulation, such that (HFUN(V2)-HFUN(V1))/LL<= 
%   DHDX, where {V1,V2} are the vertices of a given triangle
%   edge and LL is the edge-length. Limits are enforced exh-
%   austively over all edges.
%
%   See also TRIHFN2, LFSHFN2

%   This function is based on a very simplified version of: 
%   Persson, P.O. "Mesh size functions for implicit geometr-
%   ies and PDE-based gradient limiting." Engineering with 
%   Computers 22 (2006): 95-109.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 18/04/2017

%---------------------------------------------- basic checks    
    if ( ~isnumeric(vert) || ...
         ~isnumeric(tria) || ...
         ~isnumeric(hfun) || ...
         ~isnumeric(dhdx) )
        error('limhfn2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ...
        ndims(tria) ~= +2 || ...
        ndims(hfun) ~= +2 || ...
        numel(dhdx) ~= +1 )
        error('limhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || ...
        size(tria,2) < +3 || ...
        size(hfun,2)~= +1 || ...
        size(vert,1)~= size(hfun,1) )
        error('limhfn2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(vert,1) ;

%---------------------------------------------- basic checks
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nvrt )
        error('limhfn2:invalidInputArgument', ...
            'Invalid TRIA input array.') ;
    end

%-------------------- impose gradient limits over mesh edges
   [edge,tria] = tricon2(tria);
    
    evec = vert(edge(:,2),:) - ...
           vert(edge(:,1),:) ;
    elen = sqrt(sum(evec.^2,2)) ;
  
%-------------------- impose gradient limits over edge-graph  
   [hfun] = limgrad( ...
         edge,elen,hfun,dhdx,sqrt(nvrt)) ;
 
end



