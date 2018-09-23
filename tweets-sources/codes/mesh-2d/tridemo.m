function tridemo(demo)
%TRIDEMO run various triangulation demos for MESH2D.
%   TRIDEMO(N) runs the N-TH demo problem. The following de-
%   mo problems are currently available:
%
% - DEMO-0: very simple example to start with -- construct a
%   mesh for a square domain with a square hold cut from its
%   centre.
%
% - DEMO-1: explore the impact of the "radius-edge" thresho-
%   ld (RHO2) on mesh density/quality.
%
% - DEMO-2: explore the impact of the "Frontal-Delaunay" vs.
%   "Delaunay-refinement " algorithms. 
%
% - DEMO-3: explore impact of user-defined mesh-size constr-
%   aints.
%
% - DEMO-4: explore impact of "hill-climbing" mesh optimisa-
%   tions.
%
% - DEMO-5: assemble triangulations for multi-part geometry 
%   definitions.
%
% - DEMO-6: assemble triangulations for geometries with int-
%   ernal constraints.
%
% - DEMO-7: investigate the use of quadtree-type refinement.
%
% - DEMO-8: explore impact of user-defined mesh-size constr-
%   aints.
%
% - DEMO-9: larger-scale problem, mesh refinement + optimis-
%   ation. 
%
% - DEMO10: medium-scale problem, mesh refinement + optimis-
%   ation. 
%
%   See also REFINE2, SMOOTH2, TRIDIV2, FIXGEO2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 09/07/2018
%-----------------------------------------------------------

    close all;
    
    libpath();

    switch (demo)
        case  0, demo0 ();
        case  1, demo1 ();
        case  2, demo2 ();
        case  3, demo3 ();
        case  4, demo4 ();
        case  5, demo5 ();
        case  6, demo6 ();
        case  7, demo7 ();
        case  8, demo8 ();
        case  9, demo9 ();
        case 10, demo10();
            
        otherwise
        error('tridemo:invalidSelection', 'Invalid selection!') ;
    end

end

function demo0
%DEMO0 a very simple example to start with -- mesh a square
%domain with a square hold cut from its centre.

    fprintf(1, [ ...
' A very simple example to start with -- construct a mesh for \n', ...
' a simple square domain with a square hole cut from its cen- \n', ...
' tre. The geometry is specified as a Planar Straight-Line \n', ...
' Graph (PSLG) -- a list of xy coordinates, or "nodes", and a \n', ...
' list of straight-line connections between nodes, or "edges".\n', ...
' The REFINE2 routine is used to build a triangulation of the \n', ...
' domain that: (a) conforms to the geometry, and (b) contains \n', ...
' only "nicely" shaped triangles. In the second panel, a mesh \n', ...
' that additionally satisfies "mesh-size" constrains is cons- \n', ...
' structed -- '
        ] ) ;

%------------------------------------------- setup geometry
    
    node = [                % list of xy "node" coordinates
        0, 0                % outer square
        9, 0
        9, 9
        0, 9 
        4, 4                % inner square
        5, 4
        5, 5
        4, 5 ] ;
    
    edge = [                % list of "edges" between nodes
        1, 2                % outer square 
        2, 3
        3, 4
        4, 1 
        5, 6                % inner square
        6, 7
        7, 8
        8, 5 ] ;

%------------------------------------------- call mesh-gen.
   [vert,etri, ...
    tria,tnum] = refine2(node,edge) ;

%------------------------------------------- draw tria-mesh
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    
%------------------------------------------- call mesh-gen.
    hfun = +.5 ;            % uniform "target" edge-lengths

   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun) ;

%------------------------------------------- draw tria-mesh
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
end

function demo1
%DEMO1 explore impact of RHO2 threshold on mesh density/qua-
%lity.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/lake.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' The REFINE2 routine can be used to build guaranteed-quality \n', ...
' Delaunay triangulations for general polygonal geometries in \n', ...
' the two-dimensional plane. The "quality" of elements in the \n', ...
' triangulation can be controlled using the "radius-edge" bo- \n', ...
' und RHO2. \n', ...
        ] ) ;
        
%---------------------------------------------- RHO2 = +1.50   
    fprintf(1, ' \n') ;
    fprintf(1, [ ...
' Setting large values for RHO2, (RHO2 = 1.50 here) generates \n', ...
' sparse triangulations with poor worst-case angle bounds. \n', ...
        ] ) ;
    
    opts.kind = 'delaunay';
    opts.rho2 = +1.50 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['TRIA-MESH: RHO2<=+1.50, |TRIA|=' , ...
        num2str(size(tria,1))]) ;
    
%---------------------------------------------- RHO2 = +1.00    
    fprintf(1, [ ...
' Setting small values for RHO2, (RHO2 = 1.00 here) generates \n', ...
' dense triangulations with good worst-case angle bounds. \n', ...
        ] ) ;
    
    opts.kind = 'delaunay';
    opts.rho2 = +1.00 ;
    
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['TRIA-MESH: RHO2<=+1.00, |TRIA|=' , ...
        num2str(size(tria,1))]) ;
        
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
end

function demo2
%DEMO2 explore impact of refinement "KIND" on mesh quality/-
%density.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/lake.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' The REFINE2 routine supports two Delaunay-based refinement  \n', ...
' algorithms: a "standard" Delaunay-refinement type approach, \n', ...
' and a "Frontal-Delaunay" technique. For problems constrain- \n', ...
' ed by element "quality" alone, the Frontal-Delaunay approa- \n', ...
' ch typically produces sigificantly sparser meshes. in both  \n', ...
' cases, the same worst-case element quality bounds are sati- \n', ...
' fied in a guaranteed manner. \n', ...
        ] ) ;
 
%---------------------------------------------- = "DELAUNAY"
    opts.kind = 'delaunay';
    opts.rho2 = +1.00 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    hold on; axis image off;
    title(['TRIA-MESH: KIND=DELAUNAY, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
%---------------------------------------------- = "DELFRONT"
    opts.kind = 'delfront';
    opts.rho2 = +1.00 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[]  ,opts) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['TRIA-MESH: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
end

function demo3
%DEMO3 explore impact of user-defined mesh-size constraints.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/airfoil.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' Additionally, the REFINE2 routine supports size-driven ref- \n', ...
' inement, producing meshes that satisfy constraints on elem- \n', ...
' ent edge-lengths. The LFSHFN2 routine can be used to create \n', ...
' mesh-size functions based on an estimate of the "local-fea- \n', ...
' ture-size" associated with a polygonal domain. The Frontal- \n', ...
' Delaunay refinement algorithm discussed in DEMO-2 is espec- \n', ...
' ially good at generating high-quality triangulations in the \n', ...
' presence of mesh-size constraints. \n', ...
        ] ) ;
 
%---------------------------------------------- do size-fun. 
    olfs.dhdx = +0.15;
 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge, ...
                    []  ,olfs) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
    figure;
    patch('faces',tlfs(:,1:3),'vertices',vlfs , ...
        'facevertexcdata' , hlfs, ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    title(['MESH-SIZE: KIND=DELAUNAY, |TRIA|=', ...
        num2str(size(tlfs,1))]) ;
 
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['TRIA-MESH: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
end

function demo4
%DEMO4 explore impact of "hill-climbing" mesh optimisations.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/airfoil.msh'];

   [node,edge] = triread( meshfile );
 
    fprintf(1, [ ...
' The SMOOTH2 routine provides iterative mesh "smoothing" ca- \n', ...
' pabilities, seeking to improve triangulation quality by ad- \n', ...
' justing the vertex positions and mesh topology. Specifical- \n', ...
' ly, a "hill-climbing" type optimisation is implemented, gu- \n', ...
' aranteeing that mesh-quality is improved monotonically. The \n', ...
' DRAWSCR routine provides detailed analysis of triangulation \n', ...
' quality, plotting histograms of various quality metrics. \n', ...
        ] ) ;
 
%---------------------------------------------- do size-fun. 
    olfs.dhdx = +0.15;
 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge, ...
                    []  ,olfs) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-REF.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
        
%---------------------------------------------- do mesh-opt.
   [vnew,enew, ...
    tnew,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tnew(:,1:3),'vertices',vnew, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tnew,1))]) ;
           
    hvrt = trihfn2(vert,vlfs,tlfs,slfs,hlfs) ;
    hnew = trihfn2(vnew,vlfs,tlfs,slfs,hlfs) ;
    
    tricost(vert,etri,tria,tnum,hvrt) ;
    tricost(vnew,enew,tnew,tnum,hnew) ;
           
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
        
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;
    set(figure(4),'units','normalized', ...
        'position',[.35,.05,.30,.35]) ;
    
end

function demo5
%DEMO5 assemble triangulations for multi-part geometry defi-
%nitions.

    fprintf(1, [ ...
' Both the REFINE2 and SMOOTH2 routines also support "multi-  \n', ...
' part" geometry definitions -- generating conforming triang- \n', ...
' ulations that conform to internal and external constraints. \n', ...
        ] ) ;

%---------------------------------------------- create geom.
    nod1 = [
        -1., -1.; +1., -1.
        +1., +1.; -1., +1.
        ] ;
    edg1 = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
    edg1(:,3) = +0;
    
    
    nod2 = [
        +.1, +0.; +.8, +0.
        +.8, +.8; +.1, +.8
        ] ;
    edg2 = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
    edg2(:,3) = +1;
        
        
    adel = 2.*pi / +64 ;
    amin = 0.*pi ;
    amax = 2.*pi - adel;
    
    xcir = +.33 * ...
        cos(amin:adel:amax)';
    ycir = +.33 * ...
        sin(amin:adel:amax)';
    xcir = xcir - .33;
    ycir = ycir - .25;
    ncir = [xcir,ycir] ;
    numc = size(ncir,1);
        
    ecir(:,1) = ...
        [(1:numc-1)'; numc] ;
    ecir(:,2) = ...
        [(2:numc-0)'; +1  ] ;
    ecir(:,3) = +2;
    
    edg2(:,1:2) = ...
    edg2(:,1:2)+size(nod1,1);
    edge = [edg1; edg2];
    node = [nod1; nod2];
        
    ecir(:,1:2) = ...
    ecir(:,1:2)+size(node,1);
    edge = [edge; ecir];
    node = [node; ncir];
    
%-- the PART argument is a cell array that defines individu-
%-- al polygonal "parts" of the overall geometry. Each elem-
%-- ent PART{I} is a list of edge indexes, indicating which
%-- edges make up the boundary of each region.
    part{1} = [ ...
        find(edge(:,3) == 0) 
        find(edge(:,3) == 1)
        find(edge(:,3) == 2)
        ] ;
    part{2} = [ ...
        find(edge(:,3) == 1)
        ] ;
    part{3} = [ ...
        find(edge(:,3) == 2)
        ] ;
        
    edge = edge(:,1:2) ;
    
%---------------------------------------------- do size-fun.
    hmax = +0.045 ;
 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge, ...
                    part) ;
    
    hlfs = min(hmax,hlfs) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;

   [vert,etri, ...
    tria,tnum] = refine2(node,edge,part,  [], ...
                         hfun, ...
                         vlfs,tlfs,slfs,hlfs) ;
                         
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
        'facecolor',[1.,1.,1.], ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
        'facecolor',[.9,.9,.9], ...
        'edgecolor',[.2,.2,.2]) ;
    patch('faces',tria(tnum==3,1:3),'vertices',vert, ...
        'facecolor',[.8,.8,.8], ...
        'edgecolor',[.2,.2,.2]) ;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
    figure;
    patch('faces',tlfs(:,1:3),'vertices',vlfs , ...
        'facevertexcdata' , hlfs, ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    title(['MESH-SIZE: KIND=DELAUNAY, |TRIA|=', ...
        num2str(size(tlfs,1))]) ;
        
    tricost(vert,etri,tria,tnum);
           
    drawnow;
        
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ; 
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;

end

function demo6
%DEMO6 build triangulations for geometries with internal co-
%nstraints.

    fprintf(1, [ ...
' Both the REFINE2 and SMOOTH2 routines also support geometr- \n', ...
' ies containing "internal" constraints. \n', ...
        ] ) ;

%---------------------------------------------- create geom.
    node = [
        -1., -1.; +1., -1.
        +1., +1.; -1., +1.
        +.0, +.0; +.2, +.7
        +.6, +.2; +.4, +.8
        +0., +.5; -.7, +.3
        -.1, +.1; -.6, +.5
        -.9, -.8; -.6, -.7
        -.3, -.6; +.0, -.5
        +.3, -.4; -.3, +.4
        -.1, +.3
        ] ;
    edge = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
         5 ,  6 ;  5 ,  7
         5 ,  8 ;  5 ,  9
         5 , 10 ;  5 , 11
         5 , 12 ;  5 , 13
         5 , 14 ;  5 , 15
         5 , 16 ;  5 , 17
         5 , 18 ;  5 , 19
        ] ;
   
%-- the geometry must be split into its "exterior" and "int-
%-- erior" components using the optional PART argument. Each
%-- PART{I} specified should define the "exterior" boundary
%-- of a polygonal region. "Interior" constraints should not
%-- be referenced by any polygon in PART -- they are imposed
%-- as isolated edge constraints. 
    part{1} = [1,2,3,4] ;
    
%---------------------------------------------- do size-fun.
    hmax = +0.175 ;
 
%---------------------------------------------- do mesh-gen.
    opts.kind = 'delaunay' ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,part,opts, ...
                         hmax) ;
                         
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',2.0) ;
    title(['MESH-OPT.: KIND=DELAUNAY, |TRIA|=', ...
        num2str(size(tria,1))]) ;
    
    tricost(vert,etri,tria,tnum);
           
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;
        
end

function demo7
%DEMO7 investigate the use of quadtree-type mesh refinement.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/channel.msh'];

   [node,edge] = triread( meshfile );
   
    fprintf(1, [ ...
' The TRIDIV2 routine can also be used to refine existing tr- \n', ...
' angulations. Each triangle is split into four new sub-tria- \n', ...
' ngles, such that element shape is preserved. Combining the  \n', ...
' TRIDIV2 and SMOOTH2 routines allows for hierarchies of high \n', ...
' quality triangulations to be generated. \n', ...
        ] ) ;
   
%---------------------------------------------- do size-fun.
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
    pmax = max(node,[],1);
    pmin = min(node,[],1);
   
    hmax = mean(pmax-pmin)/+17 ;
    hlfs = min(hmax,hlfs);   
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
        
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
   [vnew,enew, ...
    tnew,tnum] = tridiv2(vert,etri,tria,tnum) ;
    
   [vnew,enew, ...
    tnew,tnum] = smooth2(vnew,enew,tnew,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
        
    figure;
    patch('faces',tnew(:,1:3),'vertices',vnew, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',etri(:,1:2),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tnew,1))]) ;
    
    tricost(vert,etri,tria,tnum);       
    tricost(vnew,enew,tnew,tnum);
           
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
        
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;
    set(figure(4),'units','normalized', ...
        'position',[.35,.05,.30,.35]) ;
    
end

function demo8
%DEMO8 explore impact of "hill-climbing" mesh optimisations.

%---------------------------------------------- create geom.
    node = [
        -1., -1.; +3., -1.
        +3., +1.; -1., +1.
        ] ;
    edge = [
         1 ,  2 ;  2 ,  3
         3 ,  4 ;  4 ,  1
        ] ;
        
    adel = 2.*pi / +64 ;
    amin = 0.*pi ;
    amax = 2.*pi - adel;
    
    xcir = +.20 * ...
        cos(amin:adel:amax)';
    ycir = +.20 * ...
        sin(amin:adel:amax)';
    ncir = [xcir,ycir] ;    
    numc = size(ncir,1);
        
    ecir(:,1) = ...
        [(1:numc-1)'; numc] ;
    ecir(:,2) = ...
        [(2:numc-0)'; +1  ] ;
        
    ecir = ecir+size(node,1);
    edge = [edge; ecir];
    node = [node; ncir];
    
%---------------------------------------------- do mesh-gen.
    hfun = @hfun8 ;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun);
    
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
        
    figure;
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facevertexcdata' , hfun8(vert), ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    title('MESH-SIZE function.');
   
    hvrt = feval(hfun,vert) ;
    
    tricost(vert,etri,tria,tnum,hvrt) ;
           
    drawnow;
        
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ; 
    set(figure(3),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;  
  
end

function [hfun] = hfun8(test)
%HFUN8 user-defined mesh-size function for DEMO-8.

    hmax = +.05 ;
    hmin = +.01 ;

    xmid = +0.0 ;
    ymid = +0.0 ;
    
    hcir = exp( -.5*(test(:,1)-xmid).^2 ...
                -2.*(test(:,2)-ymid).^2 ) ;

    hfun = hmax - (hmax-hmin) * hcir  ;

end

function demo9
%DEMO9 larger-scale problem, mesh refinement + optimisation.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/islands.msh'];

   [node,edge] = triread( meshfile );
 
%---------------------------------------------- do size-fun. 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
        
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
           
    tricost(vert,etri,tria,tnum);
           
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;
    
end

function demo10
%DEMO10 medium-scale problem mesh refinement + optimisation.

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    meshfile = ...
        [filepath,'/poly-data/river.msh'];

   [node,edge] = triread( meshfile );
 
%---------------------------------------------- do size-fun. 
   [vlfs,tlfs, ...
    hlfs] = lfshfn2(node,edge) ;
    
   [slfs] = idxtri2(vlfs,tlfs) ;
   
%---------------------------------------------- do mesh-gen.
    hfun = @trihfn2;
   
   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun , ...
                    vlfs,tlfs,slfs,hlfs);
        
%---------------------------------------------- do mesh-opt.
   [vert,etri, ...
    tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
           
    tricost(vert,etri,tria,tnum);
           
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.05,.05,.30,.35]) ;
    
end



