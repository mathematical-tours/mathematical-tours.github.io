function [ee,tt] = tricon2(varargin)
%TRICON2 edge-centred connectivity for a conforming 2-simpl-
%ex triangulation in the two-dimensional plane.
%   [EE,TT] = TRICON2(TT,CC) returns the edge-based adjacen-
%   cy for a mesh of 2-simlexes (triangles). EE = [V1,V2,T1,
%   T2,CE] is the set of unique 1-simplexes (edges) in the 
%   mesh TT. Each row of {V1,V2} defines an edge, each row
%   of {T1,T2} defines the two triangles adjacent to an edge 
%   and CE is a "constraint" flag, indicating which row in
%   CC (if any) the edge matches. TT = [V1,V2,V3,E1,E2,E3],
%   is the set of unique 2-simplexes in the mesh, where
%   {E1,E2,E3} define the tria-to-edge mapping. Each row of 
%   {E1,E2,E3} are the indicies of the three edges that make 
%   up each triangle.

%   Darren Engwirda : 2014 --
%   Email           : de2363@columbia.edu
%   Last updated    : 01/10/2017

%---------------------------------------------- extract args
    tt = []; cc = [];

    if (nargin>=1), tt = varargin{1}; end
    if (nargin>=2), cc = varargin{2}; end

%---------------------------------------------- basic checks
    if (~isnumeric(tt))
        error('tricon2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    if (~isnumeric(cc))
        error('tricon2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(tt) ~= +2 || size(tt,2) ~= +3)
        error('tricon2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (min(tt(:)) < +1 )
        error('tricon2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end
    
    if (~isempty(cc))
    if (ndims(cc) ~= +2 || size(cc,2) ~= +2)
        error('tricon2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    end

    isoctave = ...
    exist('OCTAVE_VERSION','builtin') > +0;
        
    nt = size(tt,1);
    nc = size(cc,1);

%------------------------------ assemble non-unique edge set
    ee = zeros(nt*3,2);
    ee((1:nt)+nt*0,:) = tt(:,[1,2]);
    ee((1:nt)+nt*1,:) = tt(:,[2,3]);
    ee((1:nt)+nt*2,:) = tt(:,[3,1]);
    
%------------------------------ unique edges and re-indexing
  %[ee, iv, jv] = ...
  %     unique(sort(ee, 2), 'rows');
   
%-- as a (much) faster alternative to the 'ROWS' based call
%-- to UNIQUE above, the edge list (i.e. pairs of UINT32 va-
%-- lues) can be cast to DOUBLE, and the sorted comparisons 
%-- performed on vector inputs! 
    ee = sort(ee,2);
   [ed,iv,jv] = unique(ee*[2^31;1]);  
    ee = ee  (iv,:);
    
%------------------- tria-to-edge indexing: 3 edges per tria 
    tt = [tt, zeros(nt*1,3)] ;
    tt(:,4) = jv((1:nt)+nt*0);
    tt(:,5) = jv((1:nt)+nt*1);
    tt(:,6) = jv((1:nt)+nt*2);
    
%------------------- edge-to-tria indexing: 2 trias per edge

    if (isoctave)
    
    %-- OCTAVE is *shockingly* bad at executing loops, so -- 
    %-- even though it involves far more operations! -- call
    %-- the vectorised version below.
    
    ne = size(ee,1);
    ee = [ee, zeros(ne*1,3)] ;
    
    ei = [tt(:,4);tt(:,5);tt(:,6)] ;
    ti = [(+1:nt),(+1:nt),(+1:nt)]';
    
   [ei,ix] = sort(ei,'ascend') ;
    ti = ti  (ix,:);
    
    ix = find(diff(ei)>=+1);
    
    ni = length(ti);
    
    ep = [+1; ix+1];
    ep = [ep; ni+1];
  
    in = ep(2:ne+1)-ep(1:ne+0) > 1 ;
  
    ee( :,3) = ti(ep(1:ne)+0);
    ee(in,4) = ti(ep(  in)+1);
    
    else
    
    %-- MATLAB is actually pretty good at JIT-ing code these
    %-- days, so use the asymptotically faster version based
    %-- on the pre-computed ordering.
    
    ne = size(ee,1);
    ee = [ee, zeros(ne*1,3)] ;
    ep = +3 * ones (ne*1,1)  ;
    for ti = +1 : nt
        ei = tt(ti,4) ; 
        ee(ei,ep(ei)) = ti;
        ej = tt(ti,5) ;
        ee(ej,ep(ej)) = ti;
        ek = tt(ti,6) ;
        ee(ek,ep(ek)) = ti;
    
        ep(ei) = ep(ei)+1 ;
        ep(ej) = ep(ej)+1 ;
        ep(ek) = ep(ek)+1 ;
    end
    
    end
    
    if (isempty(cc)), return; end
    
%------------------------------------ find constrained edges
  %[ip,ip] = ismember( ...
  %   ee(:,1:2),sort(cc,2),'rows');
   
%-- as above, the 'ROWS' based call to ISMEMBER can be sped
%-- up by casting the edge lists (i.e. pairs of UINT32 valu-
%-- es) to DOUBLE, and performing the sorted queries on vec-
%-- tor inputs!
    cc = sort(cc,2);
   [ip,ip] = ismember(ed, cc*[2^31;+1]);
   
%------------------------------------ mark constrained edges
    ee(:,5) = ip;
    
end



