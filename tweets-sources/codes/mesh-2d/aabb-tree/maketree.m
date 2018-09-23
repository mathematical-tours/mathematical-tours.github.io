function [tr] = maketree(rp,varargin)
%MAKETREE assemble an AABB search tree for a collection of 
%(hyper-)rectangles.
%   [TR] = MAKETREE(RP) returns an axis-aligned bounding-box 
%   (AABB) tree for a collection of d-rectangles embedded in 
%   R^d. The rectangles are defined as an NR-by-NDIM*2 array 
%   of min/max coordinates RP = [PMIN,PMAX], where PMIN = 
%   [A1,A2,...,ANDIM] and PMAX = [B1,B2,...,BNDIM] are the 
%   minimum/maximum coordinates associated with each rectan-
%   gle. 
%
%   The resulting tree is returned as a structure containing
%   an NT-by-NDIM*2 array of tree-node coordinates TR.XX = 
%   [PMIN,PMAX], an NT-by-2 array of tree-node indexing 
%   TR.II = [PI,CI], and an NT-by-1 cell array of node lists
%   TR.LL. PI,CI are the parent/child pointers associated 
%   with each node in the tree, and TR.LL{II} is the list of
%   rectangles associated with the II-th node.
%
%   MAKETREE forms a binary search tree, such that each 
%   "leaf" node in the tree encloses a maximum number of re-
%   ctangles. The tree is formed by recursively subdividing 
%   the bounding-box of the collection. At each division, a
%   simple heuristic is used to determine a splitting axis
%   and to position an axis-aligned splitting (hyper-)plane. 
%   The associated collection of rectangles is partitioned 
%   between two new child nodes. The dimensions of each node 
%   in the tree are selected to provide a minimal enclosure 
%   of the rectangles in its associated sub-tree. Tree nodes 
%   may overlap as a result. 
%
%   [...] = MAKETREE(RP,OP) also passes an additional user-
%   defined options structure. OP.NOBJ = {32} is the maximum 
%   allowable number of rectangles per tree-node. OP.LONG = 
%   {.75} is a relative length tolerance for "long" rectang-
%   les, such that any rectangles with RMAX(IX)-RMIN(IX) >=
%   OP.LONG * (NMAX(IX)-NMIN(IX)) remain in the parent node.
%   Here RMIN,RMAX are the coordinates of the rectangle, 
%   NMIN,NMAX are the coordinates of the enclosing node in
%   the tree, and IX is the splitting axis. Nodes that beco-
%   me "full" of "long" items may exceed their maximum rect-
%   angle capacity. OP.VTOL = {.55} is a "volume" splitting
%   criteria, designed to continue subdivision while the net
%   node volume is reducing. Specifically, a node is split 
%   if V1+V2 <= OP.VTOL*VP, where VP is the d-dim. "volume"
%   of the parent node, and V1,V2 are the volumes associated 
%   with its children.
%
%   See also DRAWTREE, QUERYSET, MAPVERT, MAPRECT

% Please see the following for additional information:
%
%   Darren Engwirda, "Locally-optimal Delaunay-refinement & 
%   optimisation-based mesh generation". Ph.D. Thesis, Scho-
%   ol of Mathematics and Statistics, Univ. of Sydney, 2014:
%   http://hdl.handle.net/2123/13148

%   Darren Engwirda : 2014 --
%   Email           : de2363@columbia.edu
%   Last updated    : 08/10/2017

    tr.xx = []; tr.ii = []; tr.ll = {}; op = [];
    
%------------------------------ quick return on empty inputs
    if (isempty(rp)), return; end
    
%---------------------------------------------- basic checks    
    if (~isnumeric(rp))
        error('maketree:incorrectInputClass', ...
            'Incorrect input class.');
    end
    
%---------------------------------------------- basic checks
    if (nargin < +1 || nargin > +2)
        error('maketree:incorrectNoOfInputs', ...
            'Incorrect number of inputs.');
    end
    if (ndims(rp) ~= +2 || ...
            mod(size(rp,2),2)~= +0)
        error('maketree:incorrectDimensions', ...
            'Incorrect input dimensions.');
    end
    
%------------------------------- extract user-defined inputs
    if (nargin >= +2), op = varargin{1}; end
    
    isoctave = exist( ...
        'OCTAVE_VERSION','builtin') > 0;
    
    if (isoctave)
    
    %-- faster for OCTAVE with large tree block size; slower
    %-- loop execution...
    
        NOBJ = +1024 ; 
        
    else
    
    %-- faster for MATLAB with small tree block size; better
    %-- loop execution... 
    
        NOBJ =  + 32 ;
        
    end
    
%--------------------------------------- user-defined inputs
    if (~isstruct(op))
        
        op.nobj = NOBJ ;
        op.long = .75;
        op.vtol = .55;
        
    else
    
    %-------------------------------- bound population count
        if (isfield(op,'nobj'))
            if (op.nobj <= +0 )
                error('Invalid options.') ;
            end
        else
            op.nobj = NOBJ ;
        end
        
    %-------------------------------- bound "long" tolerance
        if (isfield(op,'long'))
            if (op.long < +.0 || op.long > +1.)
                error('Invalid options.') ;
            end
        else
            op.long = .75; 
        end
        
    %-------------------------------- bound "long" tolerance
        if (isfield(op,'vtol'))
            if (op.vtol < +.0 || op.vtol > +1.)
                error('Invalid options.') ;
            end
        else
            op.vtol = .55; 
        end
        
    end
    
%---------------------------------- dimensions of rectangles
    nd = size(rp,2) / +2 ;
    ni = size(rp,1) ;
    
%------------------------------------------ alloc. workspace
    xl = zeros(ni*1,1*nd);
    xr = zeros(ni*1,1*nd);
    ii = zeros(ni*1,2);
    ll = cell (ni*1,1);
    ss = zeros(ni*1,1);
    
%------------------------------------ min & max coord. masks
    lv = false(size(rp,2),1);
    rv = false(size(rp,2),1);
    lv((1:nd)+nd*+0) = true ;
    rv((1:nd)+nd*+1) = true ;

%----------------------------------------- inflate rectangle
    r0 = min(rp(:,lv),[],1) ;
    r1 = max(rp(:,rv),[],1) ;
   
    rd = repmat(r1-r0,ni,1) ;
    
    rp(:,lv) = ...
    rp(:,lv) - rd * eps^.8;
    rp(:,rv) = ...
    rp(:,rv) + rd * eps^.8;
    
%----------------------------------------- rectangle centres
    rc = rp(:,lv)+rp(:,rv);
    rc = rc * .5 ;
%----------------------------------------- rectangle lengths
    rd = rp(:,rv)-rp(:,lv);

%------------------------------ root contains all rectangles
    ll{1} = (+1:ni)' ;
%------------------------------------ indexing for root node
    ii(1,1) = +0 ;
    ii(1,2) = +0 ;
%------------------------------ root contains all rectangles
    xl(1,:) = min(rp(:,lv),[],1);
    xr(1,:) = max(rp(:,rv),[],1);
    
%-- main loop : divide nodes until all constraints satisfied
    ss(+1) = +1; ns = +1; nn = +1;   
    while (ns ~= +0)
    %----------------------------------- pop node from stack
        ni = ss(ns) ; 
        ns = ns - 1 ;
    %----------------------------------- push child indexing
        n1 = nn + 1 ; 
        n2 = nn + 2 ;           
        
    %--------------------------- set of rectangles in parent
        li = ll{ni} ;    
    %--------------------------- split plane on longest axis
        dd = xr(ni,:) ...
           - xl(ni,:) ;
       [dd,ia] = sort(dd);
  
        for id = nd : -1 : +1
    %--------------------------- push rectangles to children         
            ax = ia (id) ;
            mx = dd (id) ;
    
            il = rd(li,ax) > ...
                op.long * mx ;
            lp = li( il) ;          %  "long" rectangles
            ls = li(~il) ;          % "short" rectangles
        
            if (length(lp) < ...
            0.5*length(ls)&& ...
                length(lp) < ...
            0.5 * op.nobj)
                break ;
            end 
        end
        
        if (isempty(ls) )
    %-------------------------------- partition empty, done!
            continue ;
        end
        
    % select the split position: take the mean of the set of
    % (non-"long") rectangle centres along axis AX
    %-------------------------------------------------------
        sp = sum(rc(ls,ax))/length(ls);
       
    %---------------------------- partition based on centres
        i2 = rc(ls,ax)>sp ;
        l1 = ls(~i2) ;              %  "left" rectangles
        l2 = ls( i2) ;              % "right" rectangles
    
        if (isempty(l1) || ...
            isempty(l2) )
    %-------------------------------- partition empty, done!
            continue ;
        end
        
    %-------------------------------- finalise node position
        xl(n1,:) = ...
            min(rp(l1,lv),[],1) ;
        xr(n1,:) = ...
            max(rp(l1,rv),[],1) ;
        xl(n2,:) = ...
            min(rp(l2,lv),[],1) ;
        xr(n2,:) = ...
            max(rp(l2,rv),[],1) ;
            
    %--------------------------- push child nodes onto stack        
        if (length(li) <= op.nobj )
        
            vi = prod(xr(ni,:) ...  % upper d-dim "vol."
                    - xl(ni,:) ) ;
            v1 = prod(xr(n1,:) ...  % lower d-dim "vol."
                    - xl(n1,:) ) ; 
            v2 = prod(xr(n2,:) ...
                    - xl(n2,:) ) ;
        
            if (v1+v2 < op.vtol*vi)
               
    %-------------------------------- parent--child indexing
            ii(n1,1) = ni ; 
            ii(n2,1) = ni ;
            ii(ni,2) = n1 ;
            
    %-------------------------------- finalise list indexing
            ll{ni,1} = lp ;
            ll{n1,1} = l1 ;
            ll{n2,1} = l2 ;
        
            ss(ns+1) = n1 ;
            ss(ns+2) = n2 ;
            ns = ns+2; nn = nn+2;
            
            end
            
        else
    %-------------------------------- parent--child indexing
            ii(n1,1) = ni ; 
            ii(n2,1) = ni ;
            ii(ni,2) = n1 ;
            
    %-------------------------------- finalise list indexing
            ll{ni,1} = lp ;
            ll{n1,1} = l1 ;
            ll{n2,1} = l2 ;
        
            ss(ns+1) = n1 ;
            ss(ns+2) = n2 ;
            ns = ns+2; nn = nn+2;
        
        end
        
    end
%----------------------------------------------- trim alloc.
    xl = xl(1:nn,:) ;
    xr = xr(1:nn,:) ;
    ii = ii(1:nn,:) ;
    ll = ll(1:nn,:) ;
    
%----------------------------------------------- pack struct
    tr.xx =[xl, xr] ;     
    tr.ii = ii ; 
    tr.ll = ll ;
    
end



