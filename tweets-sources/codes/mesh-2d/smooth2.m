function [vert,conn,tria,tnum] = smooth2(varargin)
%SMOOTH2 "hill-climbing" mesh-smoothing for two-dimensional,
%2-simplex triangulations.
%   [VERT,EDGE,TRIA,TNUM] = SMOOTH2(VERT,EDGE,TRIA,TNUM) re-
%   turns a "smoothed" triangulation {VERT,TRIA}, incorpora-
%   ting "optimised" vertex coordinates and mesh topology.
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
%   [VERT,EDGE,TRIA,TNUM] = SMOOTH2(... ,OPTS) passes an ad-
%   ditional options structure OPTS, containing user-defined 
%   parameters, including:
%
% - OPTS.VTOL = {+1.0E-02} -- relative vertex movement tole-
%   rance, smoothing is converged when (VNEW-VERT) <= VTOL *
%   VLEN, where VLEN is a local effective length-scale.
%
% - OPTS.ITER = {+32} -- max. number of smoothing iterations
%
% - OPTS.DISP = {+ 4} -- smoothing verbosity. Set to INF for 
%   quiet execution.
%
%   See also REFINE2, TRICOST, TRIDEMO

%   This routine is loosely based on the DISTMESH algorithm,
%   employing a "spring-based" analogy to redistribute mesh
%   vertices. Such an approach is described in: P.O. Persson 
%   and Gilbert Strang. "A simple mesh generator in MATLAB." 
%   SIAM review 46(2) 2004, pp: 329--345. Details of the al-
%   gorithm used here are somewhat different, with an alter-
%   ative spring-based update employed, in addition to hill-
%   climbing element quality guarantees, and vertex density
%   controls.

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 21/07/2017
%-----------------------------------------------------------
    
    vert = []; conn = []; tria = [] ; 
    tnum = []; 
    opts = []; hfun = []; harg = {} ;
    
%---------------------------------------------- extract args  
    if (nargin>=+1), vert = varargin{1}; end
    if (nargin>=+2), conn = varargin{2}; end
    if (nargin>=+3), tria = varargin{3}; end
    if (nargin>=+4), tnum = varargin{4}; end
    if (nargin>=+5), opts = varargin{5}; end
    
   [opts] = makeopt(opts) ;

%---------------------------------------------- default CONN
    if (isempty(conn))
      
       [edge] = tricon2(tria);

        ebnd = edge(:,4) < +1;              %-- use bnd edge
        conn = edge(ebnd,1:2);
        
    end
   
%---------------------------------------------- default TNUM
    if (isempty(tnum))
        tnum = ones(size(tria, 1), 1) ; 
    end

%---------------------------------------------- basic checks    
    if ( ~isnumeric(vert) || ...
         ~isnumeric(conn) || ...
         ~isnumeric(tria) || ...
         ~isnumeric(tnum) || ...
         ~isstruct (opts) )
        error('smooth2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ...
        ndims(conn) ~= +2 || ...
        ndims(tria) ~= +2 || ...
        ndims(tnum) ~= +2 )
        error('smooth2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (size(vert,2)~= +2 || ...
        size(conn,2)~= +2 || ...
        size(tria,2)~= +3 || ...
        size(tnum,2)~= +1 || ...
        size(tria,1)~= size(tnum,1) )
        error('smooth2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(vert,1) ;

%---------------------------------------------- basic checks
    if (min(min(conn(:,1:2))) < +1 || ...
            max(max(conn(:,1:2))) > nvrt )
        error('smooth2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nvrt )
        error('smooth2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%---------------------------------------------- output title
    if (~isinf(opts.disp))
        fprintf(1,'\n') ;
        fprintf(1,' Smooth triangulation...\n') ;
        fprintf(1,'\n') ;
        fprintf(1,[...
' -------------------------------------------------------\n', ...
'      |ITER.|          |MOVE(X)|          |DTRI(X)|     \n', ...
' -------------------------------------------------------\n', ...
             ] ) ;
    end

%---------------------------------------------- polygon bnds    
    node = vert; PSLG = conn; part = {};

    pmax = max(tnum(:));
    for ppos = +1 : pmax
    
        tsel = tnum == ppos ;
        tcur = tria(tsel,:) ;
        
       [ecur,tcur] ...
           = tricon2 (tcur) ;
    
        ebnd = ecur(:,4)==0 ;
    
        same = setset2( ...
            PSLG,ecur(ebnd,1:2));
    
        part{ppos} = find(same) ;
    
    end

%---------------------------------------------- inflate bbox
    vmin = min(vert,[],1);
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
    
%---------------------------------------------- DO MESH ITER
    tnow =  tic ;

    tcpu = struct('full',0.,'dtri',0., ...
        'tcon',0.,'iter',0.,'undo',0., ...
        'keep',0.) ;
    
    for iter = +1 : opts.iter
 
    %------------------------------------------ inflate adj.
        ttic = tic ;
       
       [edge,tria] = tricon2(tria,conn) ;
  
        tcpu.tcon = ...
            tcpu.tcon + toc(ttic) ;
    
    %------------------------------------------ compute scr.
        oscr = triscr2(vert,tria) ;
        
    %------------------------------------------ vert. iter's
        ttic =  tic ;
        
        nvrt = size(vert,1);
        nedg = size(edge,1);
        
        IMAT = sparse( ...
           edge(:,1),(1:nedg)',+1,nvrt,nedg) ;
        JMAT = sparse( ...
           edge(:,2),(1:nedg)',+1,nvrt,nedg) ;
    
        EMAT = IMAT + JMAT ;
     
        vdeg = sum(EMAT,2) ;                %-- vertex |deg|
        free = (vdeg == 0) ;
       
        vold = vert ;
        for isub = +1 : max(+2,min(+8,iter))
    
        %-- compute HFUN at vert/midpoints
            hvrt = evalhfn(vert, ...
                edge,EMAT,hfun,harg) ;
                
            hmid = hvrt(edge(:,1),:) ...
                 + hvrt(edge(:,2),:) ;
            hmid = hmid * +.5 ;
            
        %-- calc. relative edge extensions        
            evec = vert(edge(:,2),:) ...
                 - vert(edge(:,1),:) ;
            elen = ...
                sqrt(sum(evec.^2,2)) ;
            
            scal = +1.0-elen./hmid ;
            scal = min (+1.0, scal);
            scal = max (-1.0, scal);
            
        %-- projected points from each end
            ipos = vert(edge(:,1),:) ...
                 -.67*[scal,scal].*evec;
            jpos = vert(edge(:,2),:) ...
                 +.67*[scal,scal].*evec;
            
           %scal = ...                      %-- nlin. weight
           %   max(abs(scal).^.5,eps^.75); 
            scal = ...
               max(abs(scal).^ 1,eps^.75);
            
        %-- sum contributions edge-to-vert        
            vnew = ...
            IMAT*([scal,scal] .* ipos) ...
          + JMAT*([scal,scal] .* jpos) ;
      
            vsum = max(EMAT*scal,eps^.75);
            
            vnew = vnew ./ [vsum,vsum] ;

        %-- fixed points. edge projection?
            vnew(conn(:),1:2) = ...
            vert(conn(:),1:2) ;
            
            vnew(vdeg==0,1:2) = ...
            vert(vdeg==0,1:2) ;
        
        %-- reset for the next local iter.         
            vert = vnew ;
    
        end
 
        tcpu.iter = ...
            tcpu.iter + toc(ttic) ;
    
    %------------------------------------------ hill-climber
        ttic = tic ;
        
    %-- unwind vert. upadte if score lower
        nscr = ones(size(tria,1),1);
        btri = true(size(tria,1),1);
        
        umax = + 8 ;
        
        for undo = +1 : umax
        
            nscr(btri) = triscr2( ...
                vert,tria(btri,:)) ;
      
        %-- TRUE if tria needs "unwinding" 
            smin = +.70 ;
            smax = +.90 ;
            sdel = .025 ;
        
            stol = smin+iter*sdel;
            stol = min (smax,stol) ;
            
            btri = nscr <= stol ...
                 & nscr <  oscr ;
             
            if (~any(btri)), break; end
             
        %-- relax toward old vert. coord's
            ivrt = ...
               unique(tria(btri,1:3));
            
            bvrt = ...
               false(size(vert,1),1) ;
            bvrt(ivrt) = true;
            
            if (undo ~= umax)
                bnew =  +.75 ^ undo ;
                bold =  +1.0 - bnew ;
            else
                bnew =  +0.0 ;
                bold =  +1.0 - bnew ;
            end
            
            vert(bvrt,:) = ...
                bold * vold(bvrt,:) ... 
              + bnew * vert(bvrt,:) ;
      
            btri = any( ...
               bvrt(tria(:,1:3)),2) ;
        
        end
    
        oscr = nscr ;
        
        tcpu.undo = ...
            tcpu.undo + toc(ttic) ;
    
    %------------------------------------- test convergence!
        ttic = tic ;
    
        vdel = ...
            sum((vert-vold).^2,2) ;
    
        evec = vert(edge(:,2),:) ...
             - vert(edge(:,1),:) ;
        elen = ...
            sqrt(sum(evec.^2,2)) ;
        
        hvrt = evalhfn(vert, ...
            edge,EMAT,hfun,harg) ;
    
        hmid = hvrt(edge(:,1),:) ...
             + hvrt(edge(:,2),:) ;
        hmid = hmid * 0.5 ;
        scal = elen./hmid ;

        emid = vert(edge(:,1),:) ...
             + vert(edge(:,2),:) ;
        emid = emid * 0.5 ;  
      
    %------------------------------------- |deg|-based prune
        keep = false(size(vert,1),1);
        keep(vdeg>+4) = true ;
    
        keep(conn(:)) = true ;
        keep(free(:)) = true ;
        
    %------------------------------------- 'density' control
        lmax = +5. / +4.  ;
        lmin = +1. / lmax ;
        
        less = scal<=lmin ;
    	more = scal>=lmax ;
        
        vbnd = false(size(vert,1),1);
        vbnd(conn(:,1)) = true ;
        vbnd(conn(:,2)) = true ;
        
        ebad = vbnd(edge(:,1)) ...     %-- not at boundaries
             | vbnd(edge(:,2)) ;
        
        less(ebad(:)) = false;
        more(ebad(:)) = false;
              
    %------------------------------------- force as disjoint
        lidx = find (less) ;
        
        for lpos = 1 : length(lidx)
            epos = lidx(lpos,1);
            inod = edge(epos,1);
            jnod = edge(epos,2);
        %--------------------------------- if still disjoint
            if (keep(inod) && ...
                keep(jnod) )
            
            keep(inod) = false ;
            keep(jnod) = false ;

            else
                
            less(epos) = false ;
            
            end
        end
        
        ebad = ...
           keep(edge(less,1)) ...
         & keep(edge(less,2)) ;
        
        more(ebad(:)) = false ;
         
    %------------------------------------- reindex vert/tria
        redo = ...
            zeros(size(vert,1),1);
        itop = ...
            length(find(keep));
        iend = ...
            length(find(less));
        
        redo(keep) = (1:itop)';
        
        redo(edge(less,1)) = ...        %-- to new midpoints
            (itop+1 : itop+iend)';
        redo(edge(less,2)) = ...
            (itop+1 : itop+iend)';
        
        vnew =[vert(keep,:) ; 
               emid(less,:) ;
              ] ;
        tnew = redo(tria(:,1:3)) ;

        ttmp = sort(tnew,2) ;           %-- filter collapsed
        okay = all( ...
            diff(ttmp,1,2)~=0,2) ;
        okay = ...
            okay & ttmp(:,1) > 0 ;
        tnew = tnew(okay,:) ;

    %------------------------------------- quality preserver
        nscr = ...
            triscr2  (vnew,tnew) ;

        stol = +0.80 ;
        
        tbad = nscr < stol ...
             & nscr < oscr(okay) ;
    
        vbad = ...
            false(size(vnew,1),1);
        vbad(tnew(tbad,:)) = true;
             
    %------------------------------------- filter edge merge
        lidx = find (less) ;
    
        ebad = ...
        vbad(redo(edge(lidx,1))) | ...
        vbad(redo(edge(lidx,2))) ;
        
        less(lidx(ebad)) = false ;
        
        keep(edge(...
        lidx(ebad),1:2)) =  true ;
        
    %------------------------------------- reindex vert/conn
        redo = ...
            zeros(size(vert,1),1);
        itop = ...
            length(find(keep));
        iend = ...
            length(find(less));
        
        redo(keep) = (1:itop)';
        
        redo(edge(less,1)) = ...
            (itop+1 : itop+iend)';
        redo(edge(less,2)) = ...
            (itop+1 : itop+iend)';
        
        vert =[vert(keep,:);
               emid(less,:);
               emid(more,:);
              ] ;
        conn = redo(conn(:,1:2)) ;
        
        tcpu.keep = ...
            tcpu.keep + toc(ttic) ;
          
    %------------------------------------- build current CDT
        ttic = tic ;
       
       [vert,conn,tria,tnum] = ...
            deltri2 (vert, ...
                conn,node,PSLG,part) ;
          
        tcpu.dtri = ...
            tcpu.dtri + toc(ttic) ;
        
    %------------------------------------- dump-out progess!
        vdel = vdel./(hvrt.*hvrt) ;
        move = vdel > opts.vtol^2 ;
        nmov = ...
            length(find(move));
        
        ntri = size(tria,1) ;
        
        if (mod(iter,opts.disp)==+0)
            fprintf(+1, ...
            '%11i %18i %18i\n', ...
            [iter,nmov,ntri]) ;
        end
        
    %------------------------------------- loop convergence!
        if (nmov == +0), break; end
        
    end

    tria = tria( :,1:3);
 
%----------------------------------------- prune unused vert
    keep = false(size(vert,1),1);
    keep(tria(:)) = true ;
    keep(conn(:)) = true ;
  
    redo = zeros(size(vert,1),1);
    redo(keep) = ...
        (+1:length(find(keep)))';

    conn = redo(conn);
    tria = redo(tria);
    
    vert = vert(keep,:);
    
    tcpu.full = ...
        tcpu.full + toc(tnow) ;

    if (opts.dbug)
%----------------------------------------- print debug timer 
        fprintf(1,'\n') ;
        fprintf(1,' Mesh smoothing timer...\n');
        fprintf(1,'\n') ;
        fprintf(1, ...
        ' FULL: %f \n', tcpu.full);
        fprintf(1, ...
        ' DTRI: %f \n', tcpu.dtri);
        fprintf(1, ...
        ' TCON: %f \n', tcpu.tcon);
        fprintf(1, ...
        ' ITER: %f \n', tcpu.iter);
        fprintf(1, ...
        ' UNDO: %f \n', tcpu.undo);   
        fprintf(1, ...
        ' KEEP: %f \n', tcpu.keep);
        fprintf(1,'\n') ;
    end
    
    if (~isinf(opts.disp)), fprintf(1,'\n'); end

end

function [hvrt] = evalhfn(vert,edge,EMAT,hfun,harg)
%EVALHFN eval. the spacing-fun. at mesh vertices.

    if (~isempty (hfun))
        if (isnumeric(hfun))
            hvrt = hfun * ...
              ones(size(vert,1),1) ;
        else
            hvrt = feval( ...
                hfun,vert,harg{:}) ;
        end
    else
    
%-- no HFUN - HVRT is mean edge-len. at vertices!      
        evec = vert(edge(:,2),:) - ...
               vert(edge(:,1),:) ;
        elen = sqrt(sum(evec.^2,2)) ;
                
        hvrt = (EMAT*elen) ...
            ./ max(sum(EMAT,2),eps) ;
        
        free = true(size(vert,1),1) ;
        free(edge(:,1)) = false;
        free(edge(:,2)) = false;
        
        hvrt(free)      =  +inf;
        
    end

end

function [opts] = makeopt(opts)
%MAKEOPT setup the options structure for SMOOTH2.
    
    if (~isfield(opts,'iter'))
        opts.iter = +32;
    else
    if (~isnumeric(opts.iter))
        error('smooth2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.iter)~= +1)
        error('smooth2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.iter <= +0)
        error('smooth2:invalidOptionValues', ...
            'Invalid OPT.ITER selection.') ;
    end
    end
    
    if (~isfield(opts,'disp'))
        opts.disp = + 4;
    else
    if (~isnumeric(opts.disp))
        error('smooth2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.disp)~= +1)
        error('smooth2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.disp <= +0)
        error('smooth2:invalidOptionValues', ...
            'Invalid OPT.DISP selection.') ;
    end
    end
    
    if (~isfield(opts,'vtol'))
        opts.vtol = +1.0E-02;
    else
    if (~isnumeric(opts.vtol))
        error('smooth2:incorrectInputClass', ...
            'Incorrect input class.');
    end
    if (numel(opts.vtol)~= +1)
        error('smooth2:incorrectDimensions', ...
            'Incorrect input dimensions.') ;    
    end
    if (opts.vtol <= 0.)
        error('smooth2:invalidOptionValues', ...
            'Invalid OPT.VTOL selection.') ;
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



