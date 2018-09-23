function [bp,bj,tr] = findball(bb,pp,varargin)
%FINDBALL spatial queries for collections of d-balls.
%   [BP,BI] = FINDBALL(BB,PI) finds the set of d-dim. balls 
%   that intersect with a given spatial query. Balls are sp-
%   ecified as a set of centres BB(:,1:ND) and (squared)
%   radii BB(:,ND+1), where ND is the number of dimensions.
%
%   A set of intersecting balls is returned for each query 
%   point in PI, such that the II-th point is associated 
%   with the balls BI(BP(II,1):BP(II,2)). Unenclosed points 
%   have BP(II,1) == 0.
%
%   [BP,BI,TR] = FINDBALL(BB,PI) additionally returns the
%   supporting aabb-tree used internally to compute the que-
%   ry. If the underlying collection BB is static, the tree 
%   TR may be passed to subsequent calls, via [BP,BI,TR] = 
%   FINDBALL(BB,PI,TR). This syntax may lead to improved pe-
%   rformance, especially when the number of balls is large 
%   w.r.t. the number of query points. Note that in such ca-
%   ses the distribution of underlying balls is NOT permitt-
%   ed to change between calls, or erroneous results may be 
%   returned. Additional parameters used to control the cre-
%   ation of the underlying aabb-tree may also be passed via 
%   [...] = FINDBALL(BB,PI,TR,OP). See MAKETREE for additio-
%   nal information.
%
%   See also MAKETREE, FINDTRIA

% Please see the following for additional information:
%
%   Darren Engwirda, "Locally-optimal Delaunay-refinement & 
%   optimisation-based mesh generation". Ph.D. Thesis, Scho-
%   ol of Mathematics and Statistics, Univ. of Sydney, 2014:
%   http://hdl.handle.net/2123/13148

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 27/04/2017

    bp = []; bj = []; tr = []; op = [];

%---------------------------------------------- basic checks
    if (nargin < +2 || nargin > +4)
        error('findball:incorrectNumInputs', ...
            'Incorrect number of inputs.');
    end

    if (nargin >= +3), tr = varargin{1}; end
    if (nargin >= +4), op = varargin{2}; end

%------------------------------ quick return on empty inputs
    if (isempty(bb)), return; end
    
%---------------------------------------------- basic checks    
    if (~isnumeric(bb) || ~isnumeric(pp))
        error('findball:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
    if (ndims(bb) ~= +2 || size(bb,2) < +3 )
        error('findball:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    if (ndims(pp) ~= +2 || ...
            size(bb,2) ~= size(pp,2)+1)
        error('findball:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end

    if (~isempty(tr) && ~isstruct(tr) )
        error('findball:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    if (~isempty(op) && ~isstruct(op) )
        error('findball:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
      
    if (isempty(tr))
%------------------------------ compute aabb-tree for d-ball
        nd = size(pp,2) ;
        rs = sqrt(bb(:,nd+1));
        rs = rs(:,ones(1,nd));   
        ab =[bb(:,1:nd)-rs, ...     % compute aabb-tree
             bb(:,1:nd)+rs]; 
        
        tr = maketree(ab,op) ;
    end

%------------------------------ compute tree-to-vert mapping
    tm = mapvert (tr,pp);

%------------------------------ compute vert-to-ball queries
   [bi,ip,bj] = ...
        queryset (tr,tm,@ballkern,pp,bb) ;
    
%------------------------------ re-index onto full obj. list  
    bp = zeros(size(pp,1),2) ;
    bp( :,2) = -1 ;
    
    if (isempty(bi)), return ; end
    
    bp(bi,:) = ip ;

end

function [ip,ib] = ballkern(pk,bk,pp,bb)
%BALLKERN d-dim. ball-vert intersection kernel routine.

        mp = length(pk); 
        mb = length(bk);       
  
        nd = size(pp,2);

    %-------------------------- push ball/vert onto n*m tile
        pk = pk.' ;
        
        pk = pk(ones(mb,1),:); 
        pk = pk(:); 
        bk = bk(:,ones(1,mp)); 
        bk = bk(:);
  
    %-------------------------- compute O(n*m) loc. distance     
        dd = ...
    sum((pp(pk,+1:nd)-bb(bk,+1:nd)).^2,2);
        
        in = dd<=bb(bk,nd+1) ;

        ip = pk(in) ;
        ib = bk(in) ;

end



