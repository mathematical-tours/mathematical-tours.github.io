function [node,PSLG,part] = fixgeo2(varargin)
%FIXGEO2 attempts to "fix" issues with geometry definitions.
%   [NNEW,ENEW,PNEW] = FIXGEO2(NODE,EDGE,PART) returns a new
%   "repaired" geometry definition. Currently, the following
%   operations are performed:
%
%   (1) redundant nodes are "zipped" together.
%   (2) redundant edges are deleted.
%   (3) edges are split about intersecting nodes.
%   (4) edges are split about intersecting edges. 
%
%   See also REFINE2, BFSGEO2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 10/10/2017
%-----------------------------------------------------------

    node = []; PSLG = []; part = {};

%---------------------------------------------- extract ARGS
    if (nargin>=+1), node = varargin{1}; end
    if (nargin>=+2), PSLG = varargin{2}; end
    if (nargin>=+3), part = varargin{3}; end

    if (isempty(node)), return ; end

%---------------------------------------------- default EDGE
    nnum = size(node,1);
    
    if (isempty(PSLG))
        PSLG = [(1:nnum-1)',(2:nnum-0)';nnum,1];
    end
    
%---------------------------------------------- default PART    
    enum = size(PSLG,1);
    
    if (isempty(part)), part{1} = (1:enum)'; end
  
%---------------------------------------------- basic checks    
    if ( ~isnumeric(node) || ...
         ~isnumeric(PSLG) || ...
         ~iscell(part) )
        error('fixgeo2:incorrectInputClass', ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(node) ~= +2 || ...
        ndims(PSLG) ~= +2 )
        error('fixgeo2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;
    end
    if (size(node,2)~= +2 || ...
        size(PSLG,2)~= +2 )
        error('fixgeo2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;
    end
    
%---------------------------------------------- basic checks
    if (min([PSLG(:)])<+1 || ...
            max([PSLG(:)]) > nnum)
        error('fixgeo2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    pmin = cellfun(@min,part);
    pmax = cellfun(@max,part);
    
    if (min([pmin(:)])<+1 || ...
            max([pmax(:)]) > enum)
        error('fixgeo2:invalidInputs', ...
            'Invalid PART input array.') ;
    end
    
%------------------------------------- try to "fix" geometry
    while (true)

        nnum = size(node,1) ;
        enum = size(PSLG,1) ;

    %--------------------------------- prune redundant nodes
       [node,PSLG,part] = ...
            prunenode(node,PSLG,part) ;
        
    %--------------------------------- prune redundant edges
       [node,PSLG,part] = ...
            pruneedge(node,PSLG,part) ;

    %--------------------------------- node//edge intersect!
        done = false;
        while (~done)
    
       [node,PSLG,part,done] = ...
            splitnode(node,PSLG,part) ;
        
        end
        
    %--------------------------------- edge//edge intersect!
        done = false;
        while (~done)
    
       [node,PSLG,part,done] = ...
            splitedge(node,PSLG,part) ;

        end
        
        if (size(node,1) == nnum && ...
            size(PSLG,1) == enum ) 
    %--------------------------------- iterate if any change
            break ;
        end
    
    end
        
end

function [node,PSLG,part,done] ...
            = prunenode(node,PSLG,part)
%PRUNENODE "prune" redundant nodes by "zipping" those within
%tolerance of each other.

    done = true ;

%------------------------------------- calc. "zip" tolerance
    nmin = min(node,[],+1) ;
    nmax = max(node,[],+1) ;  
    ndel = nmax - nmin;
    ztol = eps  ^ 0.80;
    zlen = ztol * max(ndel);
    
%------------------------------------- index clustered nodes 
    ball = zeros(size(node,1),3) ;
    ball(:,1:2) = node(:,1:2);
    ball(:,  3) = zlen * zlen;

   [vp,vi] = ...
      findball(ball,node(:,1:2)) ;

%------------------------------------- "zip" clustered nodes
   [vt,iv] = ...
      sort(vp(:,2) - vp(:,1));
    
    izip = zeros(size(node,1),1) ;
    imap = zeros(size(node,1),1) ;
    
    for kk = size(vp,1):-1:+1
        ii = iv(kk);
        for ip = vp(ii,1) : vp(ii,2)
            jj = vi(ip) ;
            if (izip(ii) == 0 && ...
                izip(jj) == 0 && ...
                ii~= jj  )
 
            done =  false ;
            
        %----------------------------- "zip" node JJ into II
            izip(jj) = ii ;
  
            end 
        end
    end

%------------------------------------- re-index nodes//edges
    next = +1 ;
    for kk = +1:+1:size(vp,1)
        if (izip(kk) == +0)
            imap(kk) = next ;
            next = next + 1 ;
        end
    end
    
    imap(izip ~= 0) = ...
       imap(izip(izip ~= 0));
    
    PSLG = imap(PSLG) ;

    node = node(izip == 0,:);

end

function [node,PSLG,part] = pruneedge(node,PSLG,part)
%PRUNEEDGE "prune" redundant topology.

%------------------------------------- prune redundant topo. 
   [ptmp,ivec,jvec] = ...
      unique(sort(PSLG,+2),'rows') ;

    PSLG = PSLG(ivec,:);
    
    for ppos = +1 : length(part)
    
    %--------------------------------- re-index part labels!
        part{ppos} = ...
          unique(jvec(part{ppos})) ;
    
    end
    
%------------------------------------- prune collapsed topo.
    keep = diff(PSLG,[],2) ~= +0 ;
    
    jvec = zeros(size(PSLG,1),1) ;
    jvec(keep) = +1;
    jvec = cumsum(jvec);

    PSLG = PSLG(keep,:);
    
    for ppos = +1 : length(part)
    
    %--------------------------------- re-index part labels!
        part{ppos} = ...
          unique(jvec(part{ppos})) ;
    
    end

end
    
function [node,PSLG,part,done] ...
            = splitnode(node,PSLG,part)
%SPLITNODE "split" PSLG about intersecting nodes.

    done = true ;

    mark = false(size(PSLG,1),1); 
    ediv = zeros(size(PSLG,1),1);
    pair = zeros(size(PSLG,1),2);

%------------------------------------- node//edge intersect!
   [lp,li] = findline (  ...
    node(PSLG(:,1),1:2), ...
    node(PSLG(:,2),1:2),node(:,1:2)) ;
    
%------------------------------------- node//edge splitting! 
    nn = +0 ;
    for ii = +1:+1:size(lp,1)
        for ip = lp(ii,1) : lp(ii,2)
            jj = li(ip) ;
            ni = PSLG(jj,1) ;
            nj = PSLG(jj,2) ;
            if (ni~=ii && ...
                nj~=ii && ~mark(jj))

            done = false;
            
        %----------------------------- mark seen, descendent
            mark(jj) = true ;
            
            nn = nn + 1 ;
            
            pair(nn,1) = jj ;
            pair(nn,2) = ii ;
                
            end
        end
    end

    if (nn == +0), return ; end
    
%------------------------------------- re-index intersection
    pair = pair(1:nn,:);
    
    inod = PSLG(pair(:,1),1);
    jnod = PSLG(pair(:,1),2);
    
    xnod = pair(:,2);
    
    ediv(pair(:,1)) = ...
        (+1:nn)' + size(PSLG,1);
    
    PSLG(pair(:,1),1) = inod;
    PSLG(pair(:,1),2) = xnod;
    
    PSLG = [PSLG; xnod,jnod];
    
%------------------------------------- re-index edge in part
    for ppos = +1:length(part)
    
        enew = ediv(part{ppos});
        enew = enew(enew ~= 0) ;
        
        part{ppos} = ...
            [part{ppos}; enew] ;
    
    end
 
end

function [node,PSLG,part,done] ...
            = splitedge(node,PSLG,part)
%SPLITEDGE "split" PSLG about intersecting edges.

    done = true ;

    mark = zeros(size(PSLG,1),2); 
    pair = zeros(size(PSLG,1),2); 
    ediv = zeros(size(PSLG,1),1);
    flag = zeros(size(node,1),1);

%------------------------------------- edge//edge intersect!
   [lp,li] = lineline( ...
        node(PSLG(:,1),1:2), ...
        node(PSLG(:,2),1:2), ...
        node(PSLG(:,1),1:2), ...
        node(PSLG(:,2),1:2)) ;
    
%------------------------------------- edge//edge splitting!
    nn = +0 ;
%---------------------------------- parse NaN delimited data
    for ii = +1 : +1 : size(lp,1)
        
        flag(PSLG(ii,1)) = ii ;
        flag(PSLG(ii,2)) = ii ;
        
        for ip = lp(ii,1) : lp(ii,2)
            jj = li(ip) ;
            if (mark(ii) == +0 && ...
                mark(jj) == +0 && ...
                ii ~= jj)

                ni = PSLG(jj,1) ;
                nj = PSLG(jj,2) ;

                if (flag(ni) ~= ii && ...
                    flag(nj) ~= ii )

                    done = false;
                
            %------------------------- mark seen, edge-pairs
                    mark(ii)   = +1 ;
                    mark(jj)   = +1 ;
            
                    nn = nn + 1 ;

                    pair(nn,1) = ii ;
                    pair(nn,2) = jj ;

                end
                
            end
        end

    end
    
    if (nn == +0), return ; end
    
%------------------------------------- re-index intersection
    pair = pair(1:nn,:) ;
    
   [okay,tval,sval] = linenear ( ...
    node(PSLG(pair(:,1),1),1:2), ...
    node(PSLG(pair(:,1),2),1:2), ...
    node(PSLG(pair(:,2),1),1:2), ...
    node(PSLG(pair(:,2),2),1:2)) ;

    pmid = .5 * ( ...
    node(PSLG(pair(:,1),2),1:2)+ ...
    node(PSLG(pair(:,1),1),1:2)) ;
    pdel = .5 * ( ...
    node(PSLG(pair(:,1),2),1:2)- ...
    node(PSLG(pair(:,1),1),1:2)) ;

    ppos = pmid+[tval,tval].*pdel;

    qmid = .5 * ( ...
    node(PSLG(pair(:,2),2),1:2)+ ...
    node(PSLG(pair(:,2),1),1:2)) ;
    qdel = .5 * ( ...
    node(PSLG(pair(:,2),2),1:2)- ...
    node(PSLG(pair(:,2),1),1:2)) ;

    qpos = qmid+[sval,sval].*qdel;

    inod = PSLG(pair(:,1),1);
    jnod = PSLG(pair(:,1),2);
    
    anod = PSLG(pair(:,2),1);
    bnod = PSLG(pair(:,2),2);
    
    xnod = (+1:nn)'+size(node,1) ;
    
    iedg = (+1:nn)'+size(PSLG,1) ...
              + 0 * size(pair,1) ;
    jedg = (+1:nn)'+size(PSLG,1) ...
              + 1 * size(pair,1) ;
    
    ediv(pair(:,1),1) = iedg;
    ediv(pair(:,2),1) = jedg;
                      
    PSLG(pair(:,1),1) = inod;
    PSLG(pair(:,1),2) = xnod;
    
    PSLG(pair(:,2),1) = anod;
    PSLG(pair(:,2),2) = xnod;
    
    PSLG = [PSLG; xnod,jnod];
    PSLG = [PSLG; xnod,bnod];
    
    node = [node;(ppos+qpos)*.5] ;
    
%------------------------------------- re-index edge in part
    for ppos = +1:length(part)
    
        enew = ediv(part{ppos});
        enew = enew(enew ~= 0) ;
        
        part{ppos} = ...
            [part{ppos}; enew] ;
    
    end
 
end


