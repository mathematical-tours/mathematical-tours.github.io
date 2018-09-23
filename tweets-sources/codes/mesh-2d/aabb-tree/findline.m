function [lp,lj,tr] = findline(pa,pb,pp,varargin)
%FINDLINE "point-on-line" queries in d-dimensional space.
%   [LP,LI] = FINDLINE(PA,PB,PI) finds the set of d-dimensi-
%   onal line-segments that intersect with a given spatial 
%   query. Lines are specified as a set of endpoints [PA,PB]
%   where both PA and PB are NL-by-ND arrays of coordinates, 
%   where ND is the number of dimensions.
%
%   A set of intersecting lines is returned for each query 
%   point in PI, such that the II-th point is associated 
%   with the lines LI(LP(II,1):LP(II,2)). Unenclosed points 
%   have LP(II,1) == 0.
%
%   [LP,LI,TR] = FINDLINE(PA,PB,PI) additionally returns the
%   supporting aabb-tree used internally to compute the que-
%   ry. If the underlying collection [PA,PB] is static, the 
%   tree TR may be passed to subsequent calls, via 
%   [LP,LI,TR] = FINDLINE(PA,PB,PI,TR). This syntax may lead 
%   to improved performance, especially when the number of 
%   lines is large w.r.t. the number of query points. Note 
%   that in such cases the distribution of underlying lines 
%   is NOT permitted to change between calls, or erroneous 
%   results may be returned. Additional parameters used to 
%   govern the creation of the underlying aabb-tree may be 
%   passed via [...] = FINDLINE(...,TR,OP). See MAKETREE for
%   additional information.
%
%   See also MAKETREE, FINDTRIA, FINDBALL, LINELINE

% Please see the following for additional information:
%
%   Darren Engwirda, "Locally-optimal Delaunay-refinement & 
%   optimisation-based mesh generation". Ph.D. Thesis, Scho-
%   ol of Mathematics and Statistics, Univ. of Sydney, 2014:
%   http://hdl.handle.net/2123/13148

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 02/08/2017
%-----------------------------------------------------------

    lp = []; lj = []; tr = []; op = [];

%---------------------------------------------- basic checks
    if (nargin < +3 || nargin > +5)
        error('findline:incorrectNumInputs', ...
            'Incorrect number of inputs.');
    end

    if (nargin >= +4), tr = varargin{1}; end
    if (nargin >= +5), op = varargin{2}; end

%------------------------------ quick return on empty inputs
    if (isempty(pa)), return; end
    if (isempty(pb)), return; end
    
%---------------------------------------------- basic checks    
    if (~isnumeric(pa) || ~isnumeric(pb) || ...
        ~isnumeric(pp) )
        error('findline:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
    if (ndims(pa) ~= +2 || size(pa,2) < +2 || ...
        ndims(pb) ~= +2 || size(pb,2) < +2 )
        error('findline:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    if (ndims(pp) ~= +2 || ...
            size(pa,1) ~= size(pb,1) || ...
            size(pa,2) ~= size(pp,2) || ...
            size(pb,2) ~= size(pp,2) )
        error('findline:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end

    if (~isempty(tr) && ~isstruct(tr) )
        error('findline:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    if (~isempty(op) && ~isstruct(op) )
        error('findline:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
      
    if (isempty(tr))
%------------------------------ compute aabb-tree for d-line
        nd = size(pp,2) ;
        nl = size(pa,1) ;
        
        ab = zeros(nl,nd*2);  % compute aabb-tree
        for ax = +1:nd
        ab(:,ax+nd*0) = min(pa(:,ax), ...
                            pb(:,ax)) ;
        ab(:,ax+nd*1) = max(pa(:,ax), ...
                            pb(:,ax)) ;
        end        
        tr = maketree(ab,op) ;
    end

%------------------------------ compute tree-to-vert mapping
    tm = mapvert (tr,pp);

%------------------------------ compute intersection rel-tol 
    p0 = min([pa; pb],[],+1) ;
    p1 = max([pa; pb],[],+1) ;
    
    zt = max(p1-p0) * eps^.8 ;
    
%------------------------------ compute vert-to-line queries
   [li,ip,lj] = queryset( ...
       tr,tm,@linekern,pp,pa,pb,zt) ;
    
%------------------------------ re-index onto full obj. list  
    lp = zeros(size(pp,1),2) ;
    lp( :,2) = -1 ;
    
    if (isempty(li)), return ; end
    
    lp(li,:) = ip ;

end

function [ip,il] = linekern(pk,lk,pp,pa,pb,zt)
%LINEKERN d-dim. node//line intersection kernel routine.

        mp = length(pk);
        ml = length(lk);
  
    %-------------------------- push line/vert onto n*m tile
        pk = pk.' ;
        
        pk = pk(ones(ml,1),:); 
        pk = pk(:); 
        lk = lk(:,ones(1,mp)); 
        lk = lk(:);
  
    %-------------------------- compute O(n*m) intersections
        mm = (pa(lk,:)+pb(lk,:)) * +.5 ;
        DD = (pb(lk,:)-pa(lk,:)) * +.5 ;
    
        mp = mm-pp(pk,:) ;
        
        tt =-sum(mp.*DD,2)./ ...
             sum(DD.*DD,2) ;
        tt = max(min(tt,+1.),-1.);
   
        nd = size(pp,2);
   
        qq = mm + repmat(tt,1,nd) .* DD;
        
        on = ...
        sum((pp(pk,:)-qq).^2,2) <= zt^2;

        ip = pk(on) ; 
        il = lk(on) ;
        
end


