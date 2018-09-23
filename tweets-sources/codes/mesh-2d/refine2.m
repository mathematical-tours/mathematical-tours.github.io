function [vert,conn,tria,tnum] = refine2(varargin)
%REFINE2 (Frontal)-Delaunay-refinement for two-dimensional,
%polygonal geometries.
%   [VERT,EDGE,TRIA,TNUM] = REFINE2(NODE,EDGE) returns a co-
%   nstrained Delaunay triangulation of the polygonal region
%   {NODE,EDGE}. NODE is an N-by-2 array of polygonal verti-
%   ces and EDGE is an E-by-2 array of edge indexing. Each
%   row in EDGE represents an edge of the polygon, such that
%   NODE(EDGE(JJ,1),:) and NODE(EDGE(JJ,2),:) are the coord-
%   inates of the endpoints of the JJ-TH edge. If the argum-
%   ent EDGE is omitted it assumed that the vertices in NODE
%   are connected in ascending order.
%
%   [...] = REFINE2(NODE,EDGE,PART) computes a triangulation
%   for a multiply-connected geometry. PART is a cell-array 
%   of polygonal "parts", where each element PART{KK} is an 
%   array of edge indices defining a given polygonal region. 
%   EDGE(PART{KK}, :) is the set of edges in the KK-TH part.
%
%   VERT is a V-by-2 array of XY coordinates in the triangu-
%   lation, EDGE is an array of constrained edges, TRIA is a
%   T-by-3 array of triangles, and TNUM is a T-by-1 array of
%   part indices. Each row of TRIA and EDGE define an eleme-
%   nt. VERT(TRIA(II,1),:), VERT(TRIA(II,2),:) and VERT(TRIA
%   (II,3),:) are the coordinates of the II-TH triangle. The
%   edges in EDGE are defined in a similar manner. NUM is an
%   array of part indexing, such that TNUM(II) is the index 
%   of the part in which the II-TH triangle resides.
%
%   [...] = REFINE2(..., OPTS) passes an additional options 
%   structure OPTS, containing various user-defined paramet-
%   ers, including:
%
% - OPTS.KIND = {'DELFRONT'}, 'DELAUNAY' -- the type of ref-
%   inement employed. The 'DELFRONT' algorithm is typically
%   slower, but produces higher quality output.
%
% - OPTS.RHO2 = {1.025} -- the maximum allowable radius-edge 
%   ratio. Refinement proceeds until all interior triangles
%   satisfy the radius-edge threshold. Smaller radius-edge
%   ratios lead to improved triangle shape, with RHO2=1 req-
%   uiring that all angles exceed 30 degrees. Setting RHO2<1 
%   may lead to non-convergence.
%
% - OPTS.REF1 = {'REFINE'}, 'PRESERVE' -- refinement 'flag'
%   for 1-dimensional faces (i.e. edges). The 'PRESERVE' op-
%   tion results in minimal refinement, attempting to retain 
%   the initial edges without further subdivision. Edges are 
%   split only to satisfy basic geomertical conformance.
%
% - OPTS.REF2 = {'REFINE'}, 'PRESERVE' -- refinement 'flag'
%   for 2-dimensional faces (i.e. trias). The 'PRESERVE' op-
%   tion results in minimal refinement, attempting to retain 
%   the initial trias without further subdivision. Trias are 
%   split only to satisfy basic geomertical conformance. 
%
% - OPTS.SIZ1 = {1.333} -- the normalised rel.-length th-
%   reshold for edge-elements. Each exterior edge is refined 
%   until LL/HH<SIZ1, where LL is the edge-length, HH is the
%   edge-centred mesh-size value.
% 
% - OPTS.SIZ2 = {1.300} -- the normalised rel.-length th-
%   reshold for tria-elements. Each interior tria is refined
%   until RE/HH<SIZ2, where RE is an effective tria length, 
%   based on the circumradius, HH is the tria-centred mesh-
%   size value.
%
% - OPTS.DISP = { +10 } -- refinement verbosity. Set to INF
%   for quiet execution.
%
%   [...] = REFINE2(..., HFUN,ARGS) also passes an optional
%   mesh-size function argument. Setting HFUN = HMAX, where 
%   HMAX is a scalar value, imposes a constant size constra-
%   int over the full domain. HFUN can also be defined as a 
%   general function handle [HH] = HFUN(PP), where PP is an
%   N-by-2 array of XY coordinates and HH is the associated
%   vector of mesh-size values. User-defined HFUN must be
%   fully vectorised. Additional arguments {A1,A2,...AN} for 
%   HFUN can be passed as trailing parameters to REFINE2. In
%   such cases, HFUN must adopt a signature [HH] = HFUN(PP,
%   A1,A2,...,AN). HFUN must return positive values.
%
%   See also SMOOTH2, TRIDIV2, TRICOST, TRIDEMO

%   This routine implements a "multi-refinement" variant of
%   Delaunay-refinement type mesh-generation. Both standard
%   Delaunay-refinement and Frontal-Delaunay type algorithms
%   are available. The Frontal-Delaunay approach is a simpl-
%   ified version of the JIGSAW algorithm, described in:
%
% * D. Engwirda, (2014): "Locally-optimal Delaunay-refineme-
%   nt and optimisation-based mesh generation", Ph.D. Thesis 
%   School of Mathematics and Statistics, Univ. of Sydney.
%   http://hdl.handle.net/2123/13148
%
% * D. Engwirda & D. Ivers, (2016): "Off-centre Steiner poi-
%   nts for Delaunay-refinement on curved surfaces", Comput-
%   er-Aided Design, (72), 157--171.
%   http://dx.doi.org/10.1016/j.cad.2015.10.007

%   This work is an extension of the "off-centre" type tech-
%   niques introduced in: 
%
% * H. Erten & A. Ungor, (2009): "Quality triangulation with 
%   locally optimal Steiner points", SIAM Journal on Scient-
%   ific Comp. 31(3), 2103--2130.
%   http://doi.org/10.1137/080716748
%
% * S. Rebay, (1993): "Efficient Unstructured Mesh Generati-
%   on by Means of Delaunay Triangulation and Bowyer-Watson 
%   Algorithm, J. Comp. Physics 106(1), 125--138.
%   http://dx.doi.org/10.1006/jcph.1993.1097

%   Generally speaking, the Delaunay-refinement method impl-
%   emented here is a variantion of the "classical" algorit-
%   hm introduced in: 
%
% * J. Ruppert, (1995): "A Delaunay refinement algorithm for 
%   quality 2-dimensional mesh generation." Journal of Algo-
%   rithms 18(3), 548--585.
%   http://dx.doi.org/10.1006/jagm.1995.1021 
%
%   See also: S. Cheng, T. Dey & J. Shewchuk, (2012): "Dela-
%   unay mesh generation", CRC Press, for comprehensive cov-
%   erage of Delaunay-based meshing techniques.

%   A much more advanced, and fully three-dimensional imple-
%   mentation is available in the JIGSAW library. For addit-
%   ional information, see: 
%   https://github.com/dengwirda/jigsaw-matlab

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 09/07/2018
%-----------------------------------------------------------

    node = []; PSLG = []; part = {}; opts = [] ; 
    hfun = []; harg = {};

%---------------------------------------------- extract args
    if (nargin>=+1), node = varargin{1}; end
    if (nargin>=+2), PSLG = varargin{2}; end
    if (nargin>=+3), part = varargin{3}; end
    if (nargin>=+4), opts = varargin{4}; end
    if (nargin>=+5), hfun = varargin{5}; end
    if (nargin>=+6), harg = varargin(6:end); end

   [opts] = makeopt(opts) ;
 
%---------------------------------------------- default EDGE
    nnod = size(node,1) ;
    
    if (isempty(PSLG))
        PSLG = [(1:nnod-1)',(2:nnod)'; nnod,1] ;
    end
      
%---------------------------------------------- default PART    
    ncon = size(PSLG,1) ;
    
    if (isempty(part)), part{1} = (1:ncon)'; end
    
%---------------------------------------------- basic checks    
    if (~isnumeric(node) || ~isnumeric(PSLG) || ...
        ~iscell   (part) || ~isstruct (opts) )
        error('refine2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(node) ~= +2 || ndims(PSLG) ~= +2)
        error('refine2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(node,2) < +2 || size(PSLG,2) < +2)
        error('refine2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%---------------------------------------------- basic checks
    if (min([PSLG(:)])<+1 || max([PSLG(:)])>nnod)
        error('refine2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    pmin = cellfun(@min,part);
    pmax = cellfun(@max,part);
    
    if (min([pmin(:)])<+1 || max([pmax(:)])>ncon)
        error('refine2:invalidInputs', ...
            'Invalid PART input array.') ;
    end

%-------------------------------- prune any non-unique topo. 
   [ivec,ivec,jvec] = ...
        unique(sort(PSLG,+2),'rows') ;
        
    PSLG = PSLG(ivec,:) ;
    
    for ppos = +1:length(part)
    
        if ( ~isnumeric(part{ppos}) )
            error (  ...
            'refine2:incorrectInputClass', ...
                'Incorrect input class. ') ;
        end
    
        part{ppos} = ...
            unique(jvec(part{ppos})) ;
    
    end

%-------------------------------- check part "manifold-ness"
    for ppos = +1:length(part)

        eloc = PSLG(part{ppos},:) ;       
        nadj = ...
            accumarray(eloc(:),1) ;
 
        if (any(mod(nadj,2) ~= 0) )
        error('refine2:nonmanifoldInputs', ...
            'Non-manifold PART detected.') ;
        end
    
    end

%---------------------------------------------- output title
    if (~isinf(opts.disp))
        fprintf(1,'\n') ;
        fprintf(1,' Refine triangulation...\n') ;
        fprintf(1,'\n') ;
        fprintf(1,[...
' -------------------------------------------------------\n', ...
'      |ITER.|          |CDT1(X)|          |CDT2(X)|     \n', ...
' -------------------------------------------------------\n', ...
             ] ) ;
    end
    
%-------------------------------- PASS 0: inflate box bounds
    vert = node; tria = []; tnum = []; iter = 0 ;
    conn = PSLG;

    vmin = min(vert,[],1);      % inflate bbox for stability
    vmax = max(vert,[],1);
    
    vdel = vmax - 1.*vmin;
    vmin = vmin - .5*vdel;
    vmax = vmax + .5*vdel;

    vbox = [
        vmin(1), vmin(2)
        vmax(1), vmin(2)
        vmax(1), vmax(2)
        vmin(1), vmax(2)
           ] ;
    vert = [vert ; vbox] ;

%-------------------------------- PASS 0: shield sharp feat.
   [vert,conn,tria,tnum,iter] = ...
        cdtbal0(vert,conn,tria,tnum, ...
            node,PSLG,part,opts,hfun,harg,iter);

%-------------------------------- PASS 1: refine 1-simplexes
   [vert,conn,tria,tnum,iter] = ...
        cdtref1(vert,conn,tria,tnum, ...
            node,PSLG,part,opts,hfun,harg,iter);
        
%-------------------------------- PASS 2: refine 2-simplexes
   [vert,conn,tria,tnum,iter] = ...
        cdtref2(vert,conn,tria,tnum, ...
            node,PSLG,part,opts,hfun,harg,iter);
    
    if (~isinf(opts.disp)), fprintf(1,'\n'); end
        
%-------------------------------- trim extra adjacency info.
    tria = tria( :,1:3) ;
    
%-------------------------------- trim vert. - deflate bbox.
    keep = false(size(vert,1),1);
    keep(tria(:)) = true;
    keep(conn(:)) = true;
  
    redo = zeros(size(vert,1),1);
    redo(keep) = ...
        (+1:length(find(keep)))';
    
    conn = redo(conn);
    tria = redo(tria);
    
    vert = vert(keep,:) ;
    
end

function [vert,conn,tria,tnum,iter] = ...
            cdtbal0(vert,conn,tria,tnum, ...
                node,PSLG,part,opts,hfun,harg,iter)
%CDTBAL0 constrained Delaunay-refinement for "sharp" 0-dim.
%features at PSLG vertices.
%   [...] = CDTBAL0(...) refines the set of 1-simplex eleme-
%   nts incident to "sharp" features in the PSLG. Specifica-
%   lly, edges that subtend "small" angles are split about a
%   set of new "collar" vertices, equi-distributed about the
%   centre of "sharp" features. Collar size is computed as a
%   min. of the incident edge-len. and local mesh-size cons-
%   traints.

    if (iter <= opts.iter)

    %------------------------------------- build current CDT
       [vert,conn, ...
        tria,tnum] = deltri2(vert,conn, ...
                             node,PSLG, ...
                             part, ...
                             opts.dtri) ;
                            
    %------------------------------------- build current adj
       [edge,tria] = tricon2(tria,conn) ;
       
       [feat,ftri] = isfeat2(vert, ...
                             edge,tria) ;
    
        apex = false(size(vert,1), 1) ;
        apex(tria(ftri)) =  true ;
      
    %------------------------------------- eval. length-fun.
        if (~isempty(hfun))
            if (isnumeric(hfun))
            vlen = hfun * ...
              ones(size(vert,1),1) ;
            else
            vlen = feval( ...
                hfun,vert,harg{:}) ;
            vlen = vlen(:) ;
            end
        else
            vlen = +inf * ...
              ones(size(vert,1),1) ;
        end

    %------------------------------------- form edge vectors
        evec = vert(conn(:,2),:) ...
             - vert(conn(:,1),:) ;
        elen = sqrt(sum(evec.^2,2));
        evec = evec./[elen,elen] ;
        
    %------------------------------------- min. adj. lengths       
        for epos = +1 : size(conn,1)
        
            ivrt = conn(epos,1) ;
            jvrt = conn(epos,2) ;
        
            vlen(ivrt) = min( ...
            vlen(ivrt), .67*elen(epos)) ;
            vlen(jvrt) = min( ...
            vlen(jvrt), .67*elen(epos)) ;
                
        end
            
    %------------------------------------- mark feature edge     
        iref = apex(conn(:,1)) ...      %- refine at vert. 1
            & ~apex(conn(:,2)) ;
        jref = apex(conn(:,2)) ...      %- refine at vert. 2
            & ~apex(conn(:,1)) ;
        dref = apex(conn(:,1)) ...      %- refine at both!
            &  apex(conn(:,2)) ;
            
        keep =~apex(conn(:,1)) ...      %- refine at neither
            & ~apex(conn(:,2)) ;

    %------------------------------------- protecting collar
        ilen = vlen(conn(iref,1)) ;       
        inew = vert(conn(iref,1),:) ...
        + [ilen,ilen].*evec(iref,:) ;
          
        jlen = vlen(conn(jref,2)) ;       
        jnew = vert(conn(jref,2),:) ...
        - [jlen,jlen].*evec(jref,:) ;
        
        Ilen = vlen(conn(dref,1)) ;       
        Inew = vert(conn(dref,1),:) ...
        + [Ilen,Ilen].*evec(dref,:) ;
          
        Jlen = vlen(conn(dref,2)) ;       
        Jnew = vert(conn(dref,2),:) ...
        - [Jlen,Jlen].*evec(dref,:) ;
        
        vnew = [inew; jnew; Inew; Jnew] ;
     
    %------------------------------------- add new vert/edge   
        iset = (1:size(inew,1))' ...
                + size(vert,1) ;
        
        jset = (1:size(jnew,1))' ...
                + size(inew,1) + ...
                + size(vert,1) ;
                
        Iset = (1:size(Inew,1))' ...
                + size(inew,1) + ...
                + size(jnew,1) + ...
                + size(vert,1) ;
                
        Jset = (1:size(Jnew,1))' ...
                + size(inew,1) + ...
                + size(jnew,1) + ...
                + size(Inew,1) + ...
                + size(vert,1) ;
  
        vert = [vert ; vnew] ;
  
        cnew = [conn(iref,1), iset ;
                conn(iref,2), iset ;
                conn(jref,2), jset ;
                conn(jref,1), jset ;
                conn(dref,1), Iset ;
                conn(dref,2), Jset ;
                Iset, Jset] ;
        conn = [conn(keep,:); cnew ] ;
        
    end       
       
end

function [vert,conn,tria,tnum,iter] = ...
            cdtref1(vert,conn,tria,tnum, ...
                node,PSLG,part,opts,hfun,harg,iter)
%CDTREF1 constrained Delaunay-refinement for 1-simplex elem-
%nts embedded in R^2.
%   [...] = CDTREF1(...) refines the set of 1-simplex eleme-
%   nts embedded in the triangulation until all constraints 
%   are satisfied. Specifically, edges are refined until all
%   local mesh-spacing and encroachment conditions are met.
%   Refinement proceeds according to either a Delaunay-refi-
%   nement or Frontal-Delaunay type approach, depending on
%   user-settings. In either case, new steiner vertices are
%   introduced to split "bad" edges - those that violate the
%   set of prescribed constraints. In the "-DR" type process
%   edges are split about their circumballs (midpoints). In
%   the "-FD" approach, new vertices are positioned such th-
%   at mesh-spacing constraints are satisfied in a "locally-
%   optimal" fashion.

    tcpu.full = +0. ;
    tcpu.ball = +0. ;
    tcpu.hfun = +0. ;
    tcpu.encr = +0. ;
    tcpu.offc = +0. ;
    
    vidx = (1:size(vert,1))';     %- "new" vert list to test
    
    tnow =  tic ;

    ntol = +1.55;

    while (strcmpi(opts.ref1,'refine'))
        
        iter = iter + 1 ;
    
        if (iter>=opts.iter),break; end
    
    %------------------------------------- calc. circumballs
        ttic = tic ;
    
        bal1 = cdtbal1(vert,conn) ;
    
        tcpu.ball = ...
            tcpu.ball + toc(ttic) ;
        
    %------------------------------------- eval. length-fun.
        ttic = tic ;
        
        if (~isempty(hfun))
            if (isnumeric(hfun))
            fun0 = hfun * ...
              ones(size(vert,1),1);
            fun1 = hfun ;
            else
            fun0(vidx) = ...
                feval(hfun, ...
            vert(vidx,:), harg{:});
            fun0 = fun0(:) ;
            fun1 = fun0(conn(:,1))...
                 + fun0(conn(:,2));
            fun1 = fun1 / +2. ;
            end
        else
            fun0 = +inf * ...
              ones(size(vert,1),1);
            fun1 = +inf ;
        end
    
        siz1 = ...
         +4. * bal1(:,3)./(fun1.*fun1) ;
  
        tcpu.hfun = ...
            tcpu.hfun + toc(ttic) ;
  
    %------------------------------------- test encroachment
        ttic = tic ;
        
        bal1(:,3) = ...
            (1.-eps^.75) * bal1(:,3) ;
  
       [vp,vi] = ...
           findball(bal1,vert(:,1:2));

    %------------------------------------- near=>[vert,edge]
        next = +0;
        ebad = false(size(conn,1),1) ;
        near = zeros(size(conn,1),1) ;
        for ii = +1 : size(vp,1)
            for ip = vp(ii,1):vp(ii,2)
                jj = vi(ip);
                if (ii ~= conn(jj,1) ...
                &&  ii ~= conn(jj,2) )
                next = next + 1;
                near(next,1) = ii;
                near(next,2) = jj;
                end
            end
        end
        
        near = near(1:next-0,:);
        
        if (~isempty(near))
    %-- mark edge "encroached" if there is a vert within its
    %-- dia.-ball that is not joined to either of its vert's
    %-- via an existing edge...         
            ivrt = conn(near(:,2),1);
            jvrt = conn(near(:,2),2);
        
            pair = [near(:,1), ivrt];   
            ivec = setset2(pair,conn) ;
            
            pair = [near(:,1), jvrt];
            jvec = setset2(pair,conn) ;
           
            okay = ~ivec & ~jvec ;
            
            ebad(near(okay,2))=true ;
      
        end
        
        tcpu.encr = ...
            tcpu.encr + toc(ttic);
        
    %------------------------------------- refinement queues
        ref1 = false(size(conn,1),1);      
        ref1(ebad)           = true ;   %- edge encroachment
        ref1(siz1>opts.siz1* ...        %- bad equiv. length
                  opts.siz1) = true ;
        
        num1 = find(ref1)  ;
      
    %------------------------------------- dump-out progess!
        if (mod(iter,opts.disp)==0)
            numc = size(conn,1) ;
            numt = size(tria,1) ;
            fprintf(+1, ...
            '%11i %18i %18i\n', ...
            [iter,numc,numt]) ;
        end
      
    %------------------------------------- nothing to refine
        if (isempty(num1)), break; end
        
    %------------------------------------- refine "bad" tria
        switch (lower(opts.kind))
        case 'delaunay'
    %------------------------------------- do circ-ball pt's
        new1 = bal1(ref1, 1:2) ;
        
        vidx = (1:size(new1,1))' ...
                + size(vert,1) ;
        
        cnew = [conn( ref1,1), vidx
                conn( ref1,2), vidx];
        conn = [conn(~ref1,:); cnew];
        
    %------------------------------------- update vertex set    
        vert = [vert; new1(:,1:2)];
        
        
        case 'delfront'
    %-- symmetric off-centre scheme:- refine edges from both
    %-- ends simultaneously, placing new vertices to satisfy
    %-- the worst of mesh-spacing and local voronoi constra-
    %-- ints.
  
        ttic = tic ;
  
        evec = vert(conn(ref1,2),:) ...
             - vert(conn(ref1,1),:) ;
        elen = sqrt(sum(evec.^2,2)) ;
        evec = evec ./ [elen, elen] ;
  
    %------------------------------------- "voro"-type dist.
        vlen = sqrt(bal1(ref1,3));
        
    %------------------------------------- "size"-type dist.
        ihfn = fun0(conn(ref1,1));
        jhfn = fun0(conn(ref1,2));
        
    %------------------------------------- bind "safe" dist.
        ilen = min(vlen,ihfn) ;
        jlen = min(vlen,jhfn) ;
 
    %------------------------------------- locate offcentres       
        inew = vert(conn(ref1,1),:) ...
             + [ilen,ilen].*evec ;
        jnew = vert(conn(ref1,2),:) ...
             - [jlen,jlen].*evec ;
             
    %------------------------------------- iter. "size"-type
        for ioff = +1 : +3
    %------------------------------------- eval. length-fun.
        if (~isempty(hfun))
            if (isnumeric(hfun))
            iprj = hfun * ...
              ones(size(inew,1),1);
            jprj = hfun * ...
              ones(size(jnew,1),1);
            else
            iprj = feval( ...
                hfun,inew,harg{:});
            jprj = feval( ...
                hfun,jnew,harg{:});
            iprj = iprj(:);
            jprj = jprj(:);
            end
        else
            iprj = +inf * ...
              ones(size(inew,1),1);
            jprj = +inf * ...
              ones(size(jnew,1),1);
        end
        
        iprj = 0.5*ihfn + 0.5*iprj;
        jprj = 0.5*jhfn + 0.5*jprj;

    %------------------------------------- bind "safe" dist.
        ilen = min(vlen,iprj) ;
        jlen = min(vlen,jprj) ;
 
    %------------------------------------- locate offcentres       
        inew = vert(conn(ref1,1),:) ...
             + [ilen,ilen].*evec ;
        jnew = vert(conn(ref1,2),:) ...
             - [jlen,jlen].*evec ;
        
        end
        
    %------------------------------------- merge i,j if near        
        near = ...
            ilen+jlen>=vlen*ntol ;
        
        znew = inew(near,:) * .5 ...
             + jnew(near,:) * .5 ;
        
        inew = inew(~near,1:2) ;
        jnew = jnew(~near,1:2) ;
 
    %------------------------------------- split constraints
        zset = (1:size(znew,1))' ...
                + size(vert,1) ;
                
        iset = (1:size(inew,1))' ...
                + size(znew,1) + ...
                + size(vert,1) ;
                
        jset = (1:size(jnew,1))' ...
                + size(znew,1) + ...
                + size(inew,1) + ...
                + size(vert,1) ;
        
        set1 = num1( near);
        set2 = num1(~near);
        
        cnew = [conn( set1,1), zset
                conn( set1,2), zset
                conn( set2,1), iset
                conn( set2,2), jset
                iset, jset ] ;
        conn = [conn(~ref1,:); cnew];
            
    %------------------------------------- update vertex set    
        vert = [vert; znew(:,1:2)];
        vert = [vert; inew(:,1:2)];
        vert = [vert; jnew(:,1:2)];

        vidx = [zset; iset; jset] ;
        
        tcpu.offc = ...
            tcpu.offc + toc(ttic) ;
                       
           
        end % switch(lower(opts.kind))
    
    end

    tcpu.full = ...
        tcpu.full + toc(tnow) ;

    if (~isinf(opts.disp) )
    %------------------------------------- print final stats
        numc = size(conn,1) ;
        numt = size(tria,1) ;
        fprintf(+1, ...
        '%11i %18i %18i\n', ...
        [iter,numc,numt]) ;
    end

    if (opts.dbug)
    %------------------------------------- print debug timer 
        fprintf(1,'\n') ;
        fprintf(1,' 1-simplex REF. timer...\n');
        fprintf(1,'\n') ;
        fprintf(1, ...
        ' FULL: %f \n', tcpu.full);
        fprintf(1, ...
        ' BALL: %f \n', tcpu.ball);
        fprintf(1, ...
        ' HFUN: %f \n', tcpu.hfun);
        fprintf(1, ...
        ' ENCR: %f \n', tcpu.encr);
        fprintf(1, ...
        ' OFFC: %f \n', tcpu.offc);
        fprintf(1,'\n') ;
    end

end

function [vert,conn,tria,tnum,iter] = ...
            cdtref2(vert,conn,tria,tnum, ...
                node,PSLG,part,opts,hfun,harg,iter)
%CDTREF2 constrained Delaunay-refinement for 2-simplex elem-
%nts embedded in R^2.
%   [...] = CDTREF2(...) refines the set of 2-simplex eleme-
%   nts embedded in the triangulation until all constraints 
%   are satisfied. Specifically, triangles are refined until
%   all local mesh-spacing and element-shape conditions are
%   met. Refinement proceeds according to either a Delaunay-
%   refinement or Frontal-Delaunay type approach, depending 
%   on user-settings. In either case, new steiner points are
%   introduced to split "bad" triangles - those that violate 
%   the set of prescribed constraints. In the "-DR" type pr-
%   ocess triangles are split about their circumballs. In
%   the "-FD" approach, new vertices are positioned such th-
%   at mesh-spacing and element-shape constraints are satis-
%   fied in a "locally-optimal" fashion.

    tcpu.full = +0. ;
    tcpu.dtri = +0. ;
    tcpu.tcon = +0. ;
    tcpu.ball = +0. ;
    tcpu.hfun = +0. ;
    tcpu.offc = +0. ;
    tcpu.filt = +0. ;

    vidx = (1:size(vert,1))';     %- "new" vert list to test
    
    tnow =  tic ;

    near = +.775;
    
    while (strcmpi(opts.ref2,'refine'))
    
        iter = iter + 1 ;

    %------------------------------------- build current CDT
        ttic = tic ;
        
        nold = size(vert,1) ;
        
       [vert,conn, ...
        tria,tnum]= deltri2(vert,conn, ...
                            node,PSLG, ...
                            part, ....
                            opts.dtri) ;

        nnew = size(vert,1) ;
        
        vidx = ...
       [vidx; (nold:nnew)'] ;
                        
        tcpu.dtri = ...
            tcpu.dtri + toc(ttic) ;

    %------------------------------------- build current adj
        ttic = tic ;

       [edge,tria]= tricon2(tria,conn) ;
       
        tcpu.tcon = ...
            tcpu.tcon + toc(ttic) ;
        
        if (iter>=opts.iter),break; end

    %------------------------------------- calc. circumballs
        ttic = tic ;
        
        bal1 = cdtbal1(vert,conn) ;
        bal2 = cdtbal2(vert, ...
                       edge,tria) ;
        len2 = minlen2(vert,tria) ;

        rho2 = bal2(:,+3) ./ len2 ;

    %------------------------------------- refinement scores
        scr2 = rho2 .* bal2(:,+3) ;
       
        tcpu.ball = ...
            tcpu.ball + toc(ttic) ;
  
    %------------------------------------- eval. length-fun.     
        ttic = tic ;
        
        if (~isempty(hfun))
            if (isnumeric(hfun))
            fun0 = hfun * ...
              ones(size(vert,1),1);
            fun2 = hfun ;
            else
            fun0(vidx) = ...
                feval(hfun, ...
            vert(vidx,:), harg{:});
            fun0 = fun0(:) ;
            fun2 = fun0(tria(:,1))...
                 + fun0(tria(:,2))...
                 + fun0(tria(:,3));
            fun2 = fun2 / +3. ;
            end
        else
            fun0 = +inf * ...
              ones(size(vert,1),1);
            fun2 = +inf ;
        end
        
        siz2 = ...
         +3. * bal2(:,3)./(fun2.*fun2) ;

        tcpu.hfun = ...
            tcpu.hfun + toc(ttic) ;

    %------------------------------------- refinement queues
        ref1 = false(size(conn,1),1);
        ref2 = false(size(tria,1),1);
        
        stri = isfeat2(vert,edge,tria) ;
        
        ref2(rho2>opts.rho2* ...        %- bad rad-edge len.
                  opts.rho2) = true ;
        ref2(stri) = false ;
        ref2(siz2>opts.siz2* ...        %- bad equiv. length
                  opts.siz2) = true ;
        
        num2 = find(ref2);
      
    %------------------------------------- dump-out progess!
        if (mod(iter,opts.disp)==0)
            numc = size(conn,1) ;
            numt = size(tria,1) ;
            fprintf(+1, ...
            '%11i %18i %18i\n', ...
            [iter,numc,numt]) ;
        end
      
    %------------------------------------- nothing to refine
        if (isempty(num2)), break; end 
        
       [scr2,idx2] = sort( ...
            scr2(num2),'descend');
        num2 = num2(idx2);

    %------------------------------------- refine "bad" tria
        switch (lower(opts.kind))
        case 'delaunay'
    %------------------------------------- do circ-ball pt's
        new2 = zeros(length(num2),3);
        new2(:,1:2) = bal2(num2,1:2);
        
        rmin = ...                      %- min. insert radii
            len2(num2)*(1.-eps^.75)^2 ;
        
        new2(:,  3) = max( ...
            bal2(num2,3)*near^2,rmin) ;
        
        
        case 'delfront'
    %-- off-centre scheme -- refine triangles by positioning
    %-- new vertices along a local segment of the voronoi
    %-- diagram, bounded by assoc. circmballs. New points
    %-- are placed to satisfy the worst of local mesh-length 
    %-- and element-shape constraints.
     
        ttic = tic ;
    
    %------------------------------------- find frontal edge
       [lmin,emin] = ...
            minlen2(vert,tria(num2,:)) ;

        ftri = false(length(num2),1) ;
        epos = zeros(length(num2),1) ;
        tadj = zeros(length(num2),1) ;
        
        for ii = +1 : length(epos)
            epos(ii) = tria( ...
                num2(ii),emin(ii)+3) ;
        end
        
    %------------------------------------- find frontal tria
        for enum = +1 : +3
        
            eidx = tria(num2,enum+3) ;
            
            ftri = ...
            ftri | edge(eidx,5) > +0 ;
        
            ione = ...
                num2 ~= edge(eidx,3) ;
            itwo = ~ione ;
            
            tadj(ione) = ...
                edge(eidx(ione),3);
            tadj(itwo) = ...
                edge(eidx(itwo),4);
        
            okay = tadj > +0 ;
            tidx = tadj(okay);
        
            ftri(okay) = ...
            ftri(okay) | ~ref2(tidx) ;
        
        end
        
        if (~any(ftri))                 %- can this happen!?
        ftri = true(length(num2),+1) ; 
        end
       
    %------------------------------------- locate offcentres 
        emid = vert(edge(epos,+1),:) ...
             + vert(edge(epos,+2),:) ;
        emid = emid * +0.50 ;
        
        elen = sqrt(lmin(:));
        
    %------------------------------------- "voro"-type dist.    
        vvec = bal2(num2,1:2)-emid ;
        vlen = sqrt(sum(vvec.^2,2));
        vvec = vvec ./ [vlen,vlen] ;
        
        hmid = fun0(edge(epos,+1),:) ...
             + fun0(edge(epos,+2),:) ;
        hmid = hmid * +0.50 ;
        
    %------------------------------------- "ball"-type dist.
        rtri = elen * opts.off2 ;
        rfac = elen * +0.50 ;
        dsqr = rtri.^2 - rfac.^2;
        doff = rtri + ...
            sqrt(max(+0.,dsqr)) ;
        
    %------------------------------------- "size"-type dist.
        dsiz = +sqrt(3.)/2. * hmid ;
   
    %------------------------------------- bind "safe" dist.
       [dist,ioff] = ...
          min([dsiz,doff,vlen],[],2) ;
    
    %------------------------------------- locate offcentres
        off2 = ...
        emid + [dist,dist] .* vvec ;
    
    %------------------------------------- iter. "size"-type
        for isub = +1 : +3
    %------------------------------------- eval. length-fun.           
        if (~isempty(hfun))
            if (isnumeric(hfun))
            hprj = hfun * ...
              ones(size(off2,1),1) ;
            else
            hprj = feval( ...
                hfun,off2,harg{:}) ;
            hprj = hprj(:) ;
            end
        else
            hprj = +inf * ...
              ones(size(off2,1),1) ;
        end

    %------------------------------------- "size"-type dist.        
        hprj = .33*hmid + .67*hprj ;
        
        dsiz = +sqrt(3.)/2. * hprj ;
        
        dsiz(dsiz<elen*.50) = +inf ;    %- edge-ball limiter
        dsiz(dsiz>vlen*.95) = +inf ;    %- circ-ball limiter
        
    %------------------------------------- bind "safe" dist. 
       [dist,ioff] = ...
          min([dsiz,doff,vlen],[],2) ;

    %------------------------------------- locate offcentres
        off2 = ...
        emid + [dist,dist] .* vvec ;
            
        end
    
        orad = ...    
        sqrt((elen*.5).^2 + dist.^2) ;
    
    %------------------------------------- do offcentre pt's
        new2 = ...
        zeros(length(find(ftri)),+3) ;
        new2(:,1:2) = off2(ftri,1:2) ;
        
        rmin = ...                      %- min. insert radii
            lmin(ftri)*(1.-eps^.75)^2 ;
        
        new2(:,  3) = max( ...
            (orad(ftri)*near).^2,rmin);
       
        tcpu.offc = ...
            tcpu.offc + toc (ttic) ;
           
        
        end % switch(lower(opts.kind))

    %------------------------------------- inter.-ball dist.
        ttic = tic ;
      
    %------------------------------------- proximity filters
       [vp,vi] = ...
          findball(new2,new2(:,1:2)) ;
     
        keep = true (size(new2,1),1) ;      
        for ii = size(vp,1):-1:+1
            for ip = vp(ii,1) ...
                   : vp(ii,2)
                jj = vi(ip);
                if (keep(jj) && ...
                    keep(ii) && ...
                    jj < ii )
                
                keep(ii) = false ; 
                break;
          
                end 
            end
        end

        new2 = new2(keep,:);
       
    %------------------------------------- test encroachment
        bal1(:,3) = ...
            (1.-eps^.75) * bal1(:,3);
    
       [vp,vi] = ...
          findball(bal1,new2(:,1:2));
        
        keep = true (size(new2,1),1);
        for ii = +1:+1:size(vp,1)
            for ip = vp(ii,1) ...
                   : vp(ii,2)
                jj = vi(ip);
                ref1(jj) =  true ;
                keep(ii) = false ;
            end
        end
 
    %------------------------------------- leave sharp edges      
        ebnd = false(size(edge,1),1);
        ebnd(tria(stri,4:6)) = true ;
        
        enot = ...
        setset2(conn,edge(ebnd,1:2));
               
        ref1(enot) = false ;
        
    %------------------------------------- preserve boundary  
        if (strcmp(lower(opts.ref1),...
            'preserve'))
        ref1(:)    = false ;
        end
 
    %------------------------------------- refinement points
        new2 = new2(keep,:);
        new1 = bal1(ref1,:);
       
        tcpu.filt = ...
            tcpu.filt + toc(ttic) ;
        
    %------------------------------------- split constraints
        idx1 = ...
       (1:size(new1))'+size(vert,1) ;
        
        idx2 = ...
       (1:size(new2))'+size(new1,1) ...
                      +size(vert,1) ;
       
        cnew = [conn( ref1,1), idx1
                conn( ref1,2), idx1];
        conn = [conn(~ref1,:); cnew];
        
        vidx = [idx1; idx2];
 
    %------------------------------------- update vertex set              
        nold = size(vert,1);
        vert = [vert; new1(:,1:2)];
        vert = [vert; new2(:,1:2)];
        nnew = size(vert,1);

        if (nnew == nold), break; end   %- we *must* be done
        
    end

    tcpu.full = ...
        tcpu.full + toc(tnow) ;

    if (~isinf(opts.disp) )
    %------------------------------------- print final stats
        numc = size(conn,1) ;
        numt = size(tria,1) ;
        fprintf(+1, ...
        '%11i %18i %18i\n', ...
        [iter,numc,numt]) ;
    end

    if (opts.dbug)
    %------------------------------------- print debug timer 
        fprintf(1,'\n') ;
        fprintf(1,' 2-simplex REF. timer...\n');
        fprintf(1,'\n') ;
        fprintf(1, ...
        ' FULL: %f \n', tcpu.full);
        fprintf(1, ...
        ' DTRI: %f \n', tcpu.dtri);
        fprintf(1, ...
        ' TCON: %f \n', tcpu.tcon);
        fprintf(1, ...
        ' BALL: %f \n', tcpu.ball);
        fprintf(1, ...
        ' HFUN: %f \n', tcpu.hfun);
        fprintf(1, ...
        ' OFFC: %f \n', tcpu.offc);
        fprintf(1, ...
        ' FILT: %f \n', tcpu.filt);   
        fprintf(1,'\n') ;
    end

end

function [opts] = makeopt(opts)
%MAKEOPT setup the options structure for REFINE2.

    if (~isfield(opts,'dtri'))
        opts.dtri = 'constrained';
    else
    if (~strcmpi(opts.dtri, 'conforming') && ...
        ~strcmpi(opts.dtri,'constrained') )
        error( ...
    'refine2:invalidOption','Invalid constraint DTRI.'); 
    end
    end

    if (~isfield(opts,'kind'))
        opts.kind = 'delfront';
    else
    if (~strcmpi(opts.kind, 'delfront') && ...
        ~strcmpi(opts.kind, 'delaunay') )
        error( ...
    'refine2:invalidOption','Invalid refinement KIND.'); 
    end
    end
    
    if (~isfield(opts,'ref1'))
        opts.ref1 = 'refine';
    else
    if (~strcmpi(opts.ref1,   'refine') && ...
        ~strcmpi(opts.ref1, 'preserve') )
        error( ...
    'refine2:invalidOption','Invalid refinement REF1.'); 
    end
    end
    
    if (~isfield(opts,'ref2'))
        opts.ref2 = 'refine';
    else
    if (~strcmpi(opts.ref2,   'refine') && ...
        ~strcmpi(opts.ref2, 'preserve') )
        error( ...
    'refine2:invalidOption','Invalid refinement REF2.'); 
    end
    end
    
    if (~isfield(opts,'iter'))
        opts.iter = +inf;
    else
    if (~isnumeric(opts.iter))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.iter)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.iter <= +0)
        error('refine2:invalidOptionValues', ...
            'Invalid OPT.ITER selection.') ;
    end
    end
    
    if (~isfield(opts,'disp'))
        opts.disp = +10 ;
    else
    if (~isnumeric(opts.disp))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.disp)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.disp <= +0)
        error('refine2:invalidOptionValues', ...
            'Invalid OPT.DISP selection.') ;
    end
    end
    
    if (~isfield(opts,'rho2'))
        opts.rho2 = 1.025;
    else
    if (~isnumeric(opts.rho2))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.rho2)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.rho2 < +1.)
        error('refine2:invalidOptionValues', ...
            'Invalid OPT.RHO2 selection.') ;
    end
    end
    
    if (~isfield(opts,'off2'))
        opts.off2 = 0.933;
    else
    if (~isnumeric(opts.off2))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.off2)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.off2 < +.7)
        error('refine2:invalidOptionValues', ...
            'Invalid OPT.OFF2 selection.') ;
    end
    end
    
    if (~isfield(opts,'siz1'))
        opts.siz1 = 1.333;
    else
    if (~isnumeric(opts.siz1))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.siz1)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.siz1 <= 0.)
        error('refine2:invalidOptionValues', ...
            'Invalid OPT.SIZ1 selection.') ;
    end
    end
    
    if (~isfield(opts,'siz2'))
        opts.siz2 = 1.300;
    else
    if (~isnumeric(opts.siz2))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.siz2)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.siz2 <= 0.)
        error('refine2:invalidOptionValues', ...
            'Invalid OPT.SIZ2 selection.') ;
    end
    end

    if (~isfield(opts,'dbug'))
        opts.dbug = false;
    else
    if (~islogical(opts.dbug))
        error('refine2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.dbug)~= +1)
        error('refine2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    end
    
end



