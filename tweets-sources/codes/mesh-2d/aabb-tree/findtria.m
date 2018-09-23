function [tp,tj,tr] = findtria(pp,tt,pj,varargin)
%FINDTRIA spatial queries for collections of d-simplexes.
%   [TP,TI] = FINDTRIA(PP,TT,PJ) finds the set of simple-
%   xes that intersect with a given spatial query. Simplexes 
%   are specified via the vertex array PP = [X1,X2,...,XN] 
%   and the indexing array TT = [T1,T2,...,TM], such that 
%   the vertex positions for the II-th simplex are the poin-
%   ts [PP(TT(II,1),:),PP(TT(II,2),:),...,PP(TT(II,M),:)].
%
%   Simplexes are NOT required to form a conforming triangu-
%   lation. Specifically, non-delaunay, non-convex and even 
%   overlapping configurations are supported. Multiple matc-
%   hes may be returned if the collection is overlapping.
%
%   A set of intersecting simplexes is returned for each 
%   query point in PI, such that the II-th point is associa-
%   ted with the simplexes TI(TP(II,1):TP(II,2)). Unenclosed 
%   points have TP(II,1)==+0. 
%
%   In general, query points may be matched to multiple sim-
%   plexes, but in cases when single matches are guaranteed, 
%   or if only a single match is desired, the following ret-
%   urns a singly-matched indexing array that is consistent 
%   with MATLAB's existing point-location routines:
%
%      [tp,tj] = findtria(pp,tt,pj) ;
%       ti = nan(size(tp,1),1);
%       in = tp(:,1) > +0;
%       ti(in) = tj(tp(in,+1));
%
%   [TP,TI,TR] = FINDTRIA(PP,TT,PI) additionally returns the
%   supporting aabb-tree used internally to compute the que-
%   ry. If the underlying collection [PP,TT] is static, the 
%   tree TR may be passed to subsequent calls, via [...] = 
%   FINDTRIA(PP,TT,PJ,TR). This syntax may lead to improved 
%   performance, especially when the number of simplexes 
%   is large w.r.t. the number of query points. Note that in 
%   such cases the underlying simplexes are NOT permitted to 
%   change between calls, or erroneous results may be retur-
%   ned. Additional parameters used to control the creation 
%   of the underlying aabb-tree may also be passed via [...] 
%   = FINDTRIA(PP,TT,PI,TR,OP). See MAKETREE for additional 
%   information.
%
%   See also MAKETREE, QUERYSET

% Please see the following for additional information:
%
%   Darren Engwirda, "Locally-optimal Delaunay-refinement & 
%   optimisation-based mesh generation". Ph.D. Thesis, Scho-
%   ol of Mathematics and Statistics, Univ. of Sydney, 2014:
%   http://hdl.handle.net/2123/13148

%   Darren Engwirda : 2014 --
%   Email           : de2363@columbia.edu
%   Last updated    : 10/03/2018

    tp = []; tj = []; tr = []; op = [];

%---------------------------------------------- basic checks
    if (nargin < +3 || nargin > +6)
        error('queryset:incorrectNumInputs', ...
            'Incorrect number of inputs.');
    end

%------------------------------- fast return on empty inputs
    if (isempty(pj)), return; end

%------------------------------- extract user-defined inputs
    if (nargin >= +4), tr = varargin{1}; end
    if (nargin >= +5), op = varargin{2}; end
    
%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(tt) || ...
        ~isnumeric(pj) )
        error('findtria:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || size(pp,2) < +2 ...
    ||  size(pp,2) > size(tt,2))
        error('findtria:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    if (ndims(tt) ~= +2 || size(tt,2) < +3)
        error('findtria:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    
%------------------------------- test query array for VERT's
    if (ndims(pj) ~= +2 || ...
             size(pj,2) ~= size(pp,2) )
        error('findtria:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end    

%---------------------------------------------- basic checks
    if (~isempty(tr) && ~isstruct(tr) )
        error('findtria:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    if (~isempty(op) && ~isstruct(op) )
        error('findtria:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
    if (isempty(tr))
%------------------------------ compute aabb's for triangles
        bi = pp(tt(:,1),:); bj = pp(tt(:,1),:);
        for ii = 2 : size(tt,2)    
            bi = min(bi,pp(tt(:,ii),:)) ;
            bj = max(bj,pp(tt(:,ii),:)) ;
        end
        bb = [bi,bj];
    
        tr = maketree(bb,op);       % compute aabb-tree       
    end
    
%------------------------------ compute tree-to-vert mapping
    tm = mapvert (tr,pj);
    
%------------------------------ compute vert-to-tria queries
    x0 = min(pp,[],1);
    x1 = max(pp,[],1);
    rt = prod(x1 - x0) * eps^.8 ;

   [ti,ip,tj] = ...
        queryset(tr,tm,@triakern,pj,pp,tt,rt) ;
 
%------------------------------ re-index onto full obj. list  
    tp = zeros(size(pj,1),2);
    tp( :,2) = -1 ;
    
    if (isempty(ti)), return; end
    
    tp(ti,:) = ip ;
    
end

function [ip,it] = triakern(pk,tk,pi,pp,tt,rt)
%TESTPTS compute the "point-tria" matches within a tile.

    mp = length(pk); mt = length(tk);

    switch (size(tt,2))
        case 3
    %------------------------------------ pts in 2-simplexes
        pk = pk.' ;
        pk = pk(ones(mt,1),:); 
        pk = pk(:);
        tk = tk(:,ones(1,mp)); 
        tk = tk(:);

        in = intria2( ...
            pp,tt(tk,:),pi(pk,:),rt);
            
        ip = pk(in) ;
        it = tk(in) ;

        case 4
    %------------------------------------ pts in 3-simplexes
        pk = pk.' ;
        pk = pk(ones(mt,1),:); 
        pk = pk(:);
        tk = tk(:,ones(1,mp)); 
        tk = tk(:);

        in = intria3( ...
            pp,tt(tk,:),pi(pk,:),rt);
            
        ip = pk(in) ;
        it = tk(in) ;

        otherwise
    %------------------------------------ pts in d-simplexes
       [il,jl] = intrian( ...
            pp,tt(tk,:),pi(pk,:));
            
        ip = pk(il(:));
        it = tk(jl(:));
        
    end

end

function [in] = intria2(pp,tt,pi,rt)
%INTRIA2 returns TRUE for points enclosed by 2-simplexes.

    t1 = tt(:,1); t2 = tt(:,2) ; 
    t3 = tt(:,3);

    vi = pp(t1,:) - pi ;
    vj = pp(t2,:) - pi ;
    vk = pp(t3,:) - pi ;
    
%------------------------------- compute sub-volume about PI
    aa = zeros(size(tt,1),3) ;
    aa(:,1) =(vi(:,1).*vj(:,2) - ...
              vj(:,1).*vi(:,2) ) ;
    aa(:,2) =(vj(:,1).*vk(:,2) - ...
              vk(:,1).*vj(:,2) ) ;
    aa(:,3) =(vk(:,1).*vi(:,2) - ...
              vi(:,1).*vk(:,2) ) ;  
    
%------------------------------- PI is internal if same sign
    rt = rt ^ 2 ;
    in = aa(:,1).*aa(:,2) >= -rt ...
       & aa(:,2).*aa(:,3) >= -rt ...
       & aa(:,3).*aa(:,1) >= -rt ;
        
end

function [in] = intria3(pp,tt,pi,rt)
%INTRIA3 returns TRUE for points enclosed by 3-simplexes.

    t1 = tt(:,1); t2 = tt(:,2) ; 
    t3 = tt(:,3); t4 = tt(:,4) ;

    v1 = pi - pp(t1,:) ;
    v2 = pi - pp(t2,:) ;
    v3 = pi - pp(t3,:) ;
    v4 = pi - pp(t4,:) ;
    
%------------------------------- compute sub-volume about PI
    aa = zeros(size(tt,1),4) ;
    aa(:,1) = ...
   +v1(:,1).*(v2(:,2).*v3(:,3) - ...
              v2(:,3).*v3(:,2) ) ...
   -v1(:,2).*(v2(:,1).*v3(:,3) - ...
              v2(:,3).*v3(:,1) ) ...
   +v1(:,3).*(v2(:,1).*v3(:,2) - ...
              v2(:,2).*v3(:,1) ) ;
    aa(:,2) = ...
   +v1(:,1).*(v4(:,2).*v2(:,3) - ...
              v4(:,3).*v2(:,2) ) ...
   -v1(:,2).*(v4(:,1).*v2(:,3) - ...
              v4(:,3).*v2(:,1) ) ...
   +v1(:,3).*(v4(:,1).*v2(:,2) - ...
              v4(:,2).*v2(:,1) ) ;
    aa(:,3) = ...
   +v2(:,1).*(v4(:,2).*v3(:,3) - ...
              v4(:,3).*v3(:,2) ) ...
   -v2(:,2).*(v4(:,1).*v3(:,3) - ...
              v4(:,3).*v3(:,1) ) ...
   +v2(:,3).*(v4(:,1).*v3(:,2) - ...
              v4(:,2).*v3(:,1) ) ;
    aa(:,4) = ...
   +v3(:,1).*(v4(:,2).*v1(:,3) - ...
              v4(:,3).*v1(:,2) ) ...
   -v3(:,2).*(v4(:,1).*v1(:,3) - ...
              v4(:,3).*v1(:,1) ) ...
   +v3(:,3).*(v4(:,1).*v1(:,2) - ...
              v4(:,2).*v1(:,1) ) ;
              
%------------------------------- PI is internal if same sign
    rt = rt ^ 2 ;
    in = aa(:,1).*aa(:,2) >= -rt ...
       & aa(:,1).*aa(:,3) >= -rt ...
       & aa(:,1).*aa(:,4) >= -rt ...
       & aa(:,2).*aa(:,3) >= -rt ...
       & aa(:,2).*aa(:,4) >= -rt ...
       & aa(:,3).*aa(:,4) >= -rt ;
       
end

function [ii,jj] = intrian(pp,tt,pi)
%INTRIAN return a list of points and enclosing "n"-simplexes.

   [np,pd] = size(pi);
   [nt,td] = size(tt);
   
%---------------- coefficient matrices for barycentric coord.
    mm = zeros(pd,pd,nt);
    for id = +1 : pd
    for jd = +1 : pd
        mm(id,jd,:) = pp(tt(:,jd),id) - ...
                      pp(tt(:,td),id) ;
    end
    end
    
%---------------- solve linear systems for barycentric coord.
    xx = zeros(pd,np,nt);
    vp = zeros(pd,np,+1);
    for ti = +1 : nt
    %---------------------------------------- form rhs coeff.
        for id = +1 : pd
            vp(id,:) = ...
          pi(:,id) - pp(tt(ti,td),id) ;
        end
    %---------------------------- actually faster to call LU-
       [ll,uu] = lu(mm(:,:,ti));
    %---------------------------- and then forward/back solve
        xx(:,:,ti) = uu\(ll\vp);
    end  
      
%-------------------- PI is internal if coord. have same sign
    in = all(xx >= +.0-(eps^.8),1) & ...
         sum(xx,1) <= +1.+(eps^.8) ;
         
%-------------------- find lists of matching points/simplexes
   [ii,jj] = ...
        find(reshape(in, [np,nt])) ;
       
end



