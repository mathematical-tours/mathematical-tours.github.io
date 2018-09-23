function [lp,lj,tr] = lineline(pa,pb,pc,pd,varargin)
%LINELINE intersection between lines in d-dimensional space.
%   [LP,LI] = LINELINE(PA,PB,PC,PD) finds intersections bet-
%   ween line segments in d-dimensions. Lines are specified 
%   as a set of endpoints [PA,PB] and [PC,PD] where PA, PB, 
%   PC and PD are NL-by-ND arrays of coordinates, where ND 
%   is the number of dimensions.
%
%   A set of intersecting lines from [PA,PB] is returned for 
%   each query line in [PC,PB], such that the II-th query 
%   line is associated with the lines LI(LP(II,1):LP(II,2)). 
%   Lines without intersections have LP(II,1) == 0.
%
%   [LP,LI,TR] = LINELINE(PA,PB,PC,PD) additionally returns 
%   the supporting aabb-tree used internally to compute the 
%   query. If the underlying collection [PA,PB] is static, 
%   the  tree TR may be recycled for subsequent calls, using 
%   [LP,LI,TR] = FINDLINE(PA,PB,PC,PD,TR). This syntax may 
%   lead to improved performance, especially when the number 
%   of lines is large w.r.t. the number of query lines. Note 
%   that in such cases the distribution of underlying lines 
%   is NOT permitted to change between calls, or erroneous 
%   results may be returned. Additional parameters used to 
%   govern the creation of the underlying aabb-tree may be 
%   passed via [...] = LINELINE(...,TR,OP). See MAKETREE for
%   additional information.
%
%   See also MAKETREE, FINDTRIA, FINDBALL, FINDLINE

% Please see the following for additional information:
%
%   Darren Engwirda, "Locally-optimal Delaunay-refinement & 
%   optimisation-based mesh generation". Ph.D. Thesis, Scho-
%   ol of Mathematics and Statistics, Univ. of Sydney, 2014:
%   http://hdl.handle.net/2123/13148

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 10/10/2017
%-----------------------------------------------------------

    lp = []; lj = []; tr = []; op = [];

%---------------------------------------------- basic checks
    if (nargin < +4 || nargin > +6)
        error('lineline:incorrectNumInputs', ...
            'Incorrect number of inputs.');
    end

    if (nargin >= +5), tr = varargin{1}; end
    if (nargin >= +6), op = varargin{2}; end

%------------------------------ quick return on empty inputs
    if (isempty(pa)), return ; end
    if (isempty(pb)), return ; end
    if (isempty(pc)), return ; end
    if (isempty(pd)), return ; end

%---------------------------------------------- basic checks    
    if (~isnumeric(pa) || ~isnumeric(pb) || ...
        ~isnumeric(pc) || ~isnumeric(pd) )
        error('lineline:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
    if (ndims(pa) ~= +2 || size(pa,2) < +2 || ...
        ndims(pb) ~= +2 || size(pb,2) < +2 || ...
        ndims(pc) ~= +2 || size(pc,2) < +2 || ...
        ndims(pd) ~= +2 || size(pd,2) < +2 )
        error('lineline:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    if (size(pa,1) ~= size(pb,1) || ...
        size(pc,1) ~= size(pd,1) )
        error('lineline:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end

    if (~isempty(tr) && ~isstruct(tr) )
        error('lineline:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    if (~isempty(op) && ~isstruct(op) )
        error('lineline:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
    nd = size(pa,2);
    nl = size(pa,1);
    ml = size(pc,1);
    
    if (isempty(tr))
%------------------------------ compute aabb-tree for d-line
    ab = zeros(nl,nd*+2) ;
    for ax = +1:nd            % compute aabb's
    ab(:,ax+nd*0) = min(pa(:,ax), ...
                        pb(:,ax)) ;
    ab(:,ax+nd*1) = max(pa(:,ax), ...
                        pb(:,ax)) ;
    end        
    tr = maketree(ab,op) ;  
    
    end

%------------------------------ compute tree-to-vert mapping
    ab = zeros(ml,nd*+2) ;
    for ax = +1:nd            % compute aabb's
    ab(:,ax+nd*0) = min(pc(:,ax), ...
                        pd(:,ax)) ;
    ab(:,ax+nd*1) = max(pc(:,ax), ...
                        pd(:,ax)) ;
    end        
    tm = maprect (tr,ab) ;
    
%------------------------------ compute line-to-line queries
   [li,ip,lj] = queryset( ...
    tr,tm,@linekern,pc,pd,pa,pb) ;
    
%------------------------------ re-index onto full obj. list  
    lp = zeros(size(pc,1),2) ;
    lp( :,2) = -1 ;
    
    if (isempty(li)), return ; end
    
    lp(li,:) = ip ;
    
end

function [i1,i2] = linekern(l1,l2,pa,pb,pc,pd)
%LINEKERN d-dim. line//line intersection kernel routine.

        m1 = length(l1) ;
        m2 = length(l2) ;
  
    %-------------------------- push line/vert onto n*m tile
        l1 = l1.' ;
        
        l1 = l1(ones(m2,1),:) ; 
        l1 = l1(:); 
        l2 = l2(:,ones(1,m1)) ; 
        l2 = l2(:);
  
    %-------------------------- compute O(n*m) intersections
       [ok,tp,tq] = linenear( ...
           pa(l1,:),pb(l1,:), ...
           pc(l2,:),pd(l2,:)) ;
    
        rt = +1.+eps;
       
        ix = abs(tp) <= +rt & ...
             abs(tq) <= +rt ;
       
        ix(~ok) = false ;
         
        i1 = l1(ix) ; 
        i2 = l2(ix) ;
    
end


