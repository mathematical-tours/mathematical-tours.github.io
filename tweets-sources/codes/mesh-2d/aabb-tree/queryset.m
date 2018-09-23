function [qi,qp,qj] = queryset(tr,tm,fn,varargin)
%QUERYSET spatial queries for AABB-indexed collections. 
%   [QI,QP,QJ] = QUERYSET(TR,TM,FN) computes a spatial query
%   on an indexed collection. TR is the AABB-tree built to 
%   index the collection, TM is the query-to-tree mapping
%   structure, and FN is the intersection "kernel" function,
%   called to compute actual intersections.
%
%   A set of intersecting objects is returned for each item, 
%   such that for each query item QI(II), a list of interse-
%   cting objects QJ(QP(II,1):QP(II,2)) is returned. 
%
%   [PK,CK] = FN(PJ,CJ,A1,...,AN) is an intersection kernel 
%   function called for each non-empty node in the tree TR. 
%   PJ,CJ are Nx1 and Mx1 arrays of query- and obj.-indices 
%   to compare against each other. These lists represent a 
%   "localised" O(N*M) "tile" of pairwise comparisons. PK,CK
%   are lists of matching objects. A1,...,AN are a set of 
%   optional user-defined arguments that are passed to the 
%   kernel function FN internally. Optional arguments can be 
%   passed to FN via calls to QUERYSET(..., A1,...AN).
%
%   See also MAPVERT, MAPRECT, MAKETREE

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 07/10/2017

    qi = []; qp = []; qj = [];

%---------------------------------------------- basic checks
    if (nargin <= +2)
        error('queryset:incorrectNumInputs', ...
            'Incorrect number of inputs.');
    end

%---------------------------------------------- empty inputs
    if (isempty(tr)), return; end

%---------------------------------------------- basic checks    
    if (~isempty(tr) && ~isstruct(tr) )
        error('queryset:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    if (~isempty(tm) && ~isstruct(tm) )
        error('queryset:incorrectInputClass', ...
            'Incorrect input class.') ;
    end

%---------------------------------- check existing aabb-tree
    if (~isfield(tm,'ii') || ...
        ~isfield(tm,'ll') )
        error('queryset:incorrectAABBstruct', ...
            'Invalid aabb-maps obj.') ;
    end
%---------------------------------- check existing aabb-tree
    if (~isfield(tr,'xx') || ...
        ~isfield(tr,'ii') || ...
        ~isfield(tr,'ll') )
        error('queryset:incorrectAABBstruct', ...
            'Invalid aabb-tree obj.') ;
    end
    
%------------------------------ spatial query over tree-node
    ic = cell(size(tm.ii,1),1) ;
    jc = cell(size(tm.ii,1),1) ;
    
    for ip = 1 : size(tm.ii,1)
    %-------------------------- extract items/query per tile
        ni = tm.ii(ip,1) ;          % node (in tree)
        
    %-------------------------- do O(n*m) search within tile
       [qi,qj] = feval(fn, ...
            tm.ll{ip,1}, ...        % query in tile
            tr.ll{ni,1}, ...        % items in tile
            varargin {:} ) ;
        
    %-------------------------- push loc. item-query matches
        ic{ip} = qi(:) ;
        jc{ip} = qj(:) ;
    end
    
%-------------------------------- concat matches into arrays
    qi = vertcat(ic{:});
    qj = vertcat(jc{:});
    
    if (isempty(qj)),return; end
    
%-------------------------------- form sparse-style indexing
   [qi,ix] = sort (qi) ; 
    qj = qj(ix);
    ix = find(diff(qi));
    
    ni = length (qi) ;
    
    qi = qi([ix;ni]) ;
    
    nj = length (qj) ;
    ni = length (qi) ;
    
%------------------------------ each list is IP(I,1):IP(I,2)
    qp = zeros(ni,2) ;
    qp(:,1) = [+1; ix+1] ;
    qp(:,2) = [ix; nj+0] ;

end



