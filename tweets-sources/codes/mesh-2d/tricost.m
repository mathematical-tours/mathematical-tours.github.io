function tricost(varargin)
%TRICOST draw quality-metrics for a 2-simplex triangulation
%embedded in the two-dimensional plane.
%   TRICOST(VERT,EDGE,TRIA,TNUM) draws histograms of quality
%   metrics for the triangulation.
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
%   TRICOST(...,HVRT) additionally draws histograms of rela-
%   tive edge-length, indicating conformance to the spacing
%   constraints. HVRT is a V-by-1 array of spacing informat-
%   ion, per an evaluation of the mesh-size function at the
%   mesh vertices VERT.
%
%   See also REFINE2, SMOOTH2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 09/07/2018
%-----------------------------------------------------------

    vert = [] ; conn = [] ; tria = [] ; 
    tnum = [] ; hvrt = [] ;

%---------------------------------------------- extract args
    if (nargin>=+1), vert = varargin{1}; end
    if (nargin>=+2), conn = varargin{2}; end
    if (nargin>=+3), tria = varargin{3}; end
    if (nargin>=+4), tnum = varargin{4}; end
    if (nargin>=+5), hvrt = varargin{5}; end

%---------------------------------------------- basic checks    
    if ( ~isnumeric(vert) || ...
         ~isnumeric(conn) || ...
         ~isnumeric(tria) || ...
         ~isnumeric(tnum) || ...
         ~isnumeric(hvrt) )
        error('tricost:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ...
        ndims(conn) ~= +2 || ...
        ndims(tria) ~= +2 )
        error('tricost:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || ...
        size(conn,2) < +2 || ...
        size(tria,2) < +3 )
        error('tricost:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(vert,1) ;
    ntri = size(tria,1) ;

%---------------------------------------------- basic checks
    if (min(min(conn(:,1:2))) < +1 || ...
            max(max(conn(:,1:2))) > nvrt )
        error('tricost:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
 
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nvrt )
        error('tricost:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%-- borrowed from the JIGSAW library!

%-- draw sub-axes directly -- sub-plot gives
%-- silly inconsistent ax spacing...!
    
    axpos31 = [.125,.750,.800,.150] ;
    axpos32 = [.125,.450,.800,.150] ;
    axpos33 = [.125,.150,.800,.150] ;
    
    axpos41 = [.125,.835,.800,.135] ;
    axpos42 = [.125,.590,.800,.135] ;
    axpos43 = [.125,.345,.800,.135] ;
    axpos44 = [.125,.100,.800,.135] ;
    
%-- draw cost histograms for 2-tria elements
    figure;
    set(gcf,'color','w','units','normalized', ...
        'position',[.05,.10,.30,.30]);
    if (~isempty(hvrt))
    
%-- have size-func data
    axes('position',axpos41); hold on;
    scrhist(triscr2(vert,tria),'tria3');
    axes('position',axpos42); hold on;
    anghist(triang2(vert,tria),'tria3');
    axes('position',axpos43); hold on;
    hfnhist(relhfn2(vert, ...
                    tria,hvrt),'tria3');
    axes('position',axpos44); hold on;
    deghist(trideg2(vert,tria),'tria3');
    
    else
    
%-- null size-func data
    axes('position',axpos31); hold on;
    scrhist(triscr2(vert,tria),'tria3');
    axes('position',axpos32); hold on;
    anghist(triang2(vert,tria),'tria3');
    axes('position',axpos33); hold on;
    deghist(trideg2(vert,tria),'tria3');
    
    end
      
end

function [mf] = mad(ff)
%MAD return mean absolute deviation (from the mean).

    mf = mean(abs(ff-mean(ff))) ;
    
end

function deghist(dd,ty)
%DEGHIST draw histogram for "degree" quality-metric.

    dd = dd(:);
    be = 1:max(dd);
    hc = histc(dd,be);
    
    r = [.85,.00,.00] ; y = [1.0,.95,.00] ;
    g = [.00,.90,.00] ; k = [.60,.60,.60] ;
    
    bar(be,hc,1.05,'facecolor',k,'edgecolor',k);
    
    axis tight;
    set(gca,'ycolor', get(gca,'color'),'ytick',[],...
        'xtick',2:2:12,'layer','top','fontsize',...
            14,'linewidth',2.,'ticklength',[.025,.025],...
                'box','off','xlim',[0,12]);
     
    switch (ty)
    case 'tria4'
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-.225,0,'$|d|_{\tau}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-.225,0, '|d|_{\tau}' ,...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
        
    case 'tria3'
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-.225,0,'$|d|_{f}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-.225,0, '|d|_{\tau}' ,...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
        
    end
    
end

function anghist(ad,ty)
%ANGHIST draw histogram for "angle" quality-metric.

    ad = ad(:);
    be = linspace(0.,180.,91);
    bm =(be(1:end-1)+be(2:end))/2.;
    hc = histc(ad,be);
    
    switch (ty)
    case 'tria4'
        poor = bm <  10.  | bm >= 160. ;
        okay =(bm >= 10.  & bm <  20. )| ...
              (bm >= 140. & bm <  160.);
        good =(bm >= 20.  & bm <  30. )| ...
              (bm >= 120. & bm <  140.);
        best = bm >= 30.  & bm <  120. ;
        
    case 'tria3'
        poor = bm <  15.  | bm >= 150. ;
        okay =(bm >= 15.  & bm <  30. )| ...
              (bm >= 120. & bm <  150.);
        good =(bm >= 30.  & bm <  45. )| ...
              (bm >= 90.  & bm <  120.);
        best = bm >= 45.  & bm <  90.  ;
    
    end
    
    r = [.85,.00,.00] ; y = [1.0,.95,.00] ;
    g = [.00,.90,.00] ; k = [.60,.60,.60] ;
    
    bar(bm(poor),hc(poor),1.05,...
        'facecolor',r,'edgecolor',r) ;
    bar(bm(okay),hc(okay),1.05,...
        'facecolor',y,'edgecolor',y) ;
    bar(bm(good),hc(good),1.05,...
        'facecolor',g,'edgecolor',g) ;
    bar(bm(best),hc(best),1.05,...
        'facecolor',k,'edgecolor',k) ;
    
    axis tight;
    set(gca,'ycolor', get(gca,'color'),'ytick',[],...
        'xtick',0:30:180,'layer','top','fontsize',...
            14,'linewidth',2.,'ticklength',[.025,.025],...
                'box','off','xlim',[0.,180.]) ;
    
    mina = max(1.000,min(ad)); %%!! so that axes don't obscure!
    maxa = min(179.0,max(ad));
    
    bara = mean(ad(:));
    mada = mad (ad(:));
    
    line([ mina, mina],...
        [0,max(hc)],'color','r','linewidth',1.5);
    line([ maxa, maxa],...
        [0,max(hc)],'color','r','linewidth',1.5);
    
    if ( mina > 25.0)
        text(mina-1.8,.90*max(hc),num2str(min(ad),'%16.1f'),...
            'horizontalalignment',...
                'right','fontsize',15) ;
    else
        text(mina+1.8,.90*max(hc),num2str(min(ad),'%16.1f'),...
            'horizontalalignment',...
                'left' ,'fontsize',15) ;
    end
    
    if ( maxa < 140.)
        text(maxa+1.8,.90*max(hc),num2str(max(ad),'%16.1f'),...
            'horizontalalignment',...
                'left' ,'fontsize',15) ; 
    else
        text(maxa-1.8,.90*max(hc),num2str(max(ad),'%16.1f'),...
            'horizontalalignment',...
                'right','fontsize',15) ;     
    end
    
    if ( maxa < 100.)
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(maxa-16.,.45*max(hc),...
        '$\bar{\sigma}_{\theta}\!= $',...
            'horizontalalignment', 'left',...
                'fontsize',16,'interpreter','latex') ;
    
        text(maxa+1.8,.45*max(hc),num2str(mad(ad),'%16.2f'),...
            'horizontalalignment',...
                'left' ,'fontsize',15) ;             
    end
                
    else
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(maxa-16.,.45*max(hc),...
        '$\bar{\sigma}_{\theta}\!= $',...
            'horizontalalignment', 'left',...
                'fontsize',16,'interpreter','latex') ;
    
        text(maxa+1.8,.45*max(hc),num2str(mad(ad),'%16.3f'),...
            'horizontalalignment',...
                'left' ,'fontsize',15) ;
    end
                
    end
    
    switch (ty)
    case 'tria4'
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-9.0,0.0,'$\theta_{\tau}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-9.0,0.0, '\theta_{\tau}' ,...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
        
    case 'tria3'
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-9.0,0.0,'$\theta_{f}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-9.0,0.0, '\theta_{f}' ,...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
        
    end
    
end

function scrhist(sc,ty)
%SCRHIST draw histogram for "score" quality-metric.

    sc = sc(:);
    be = linspace(0.,1.,101);
    bm = (be(1:end-1)+be(2:end)) / 2.;
    hc = histc(sc,be);

    switch (ty)   
    case{'tria4','dual4'}
        poor = bm <  .25 ;
        okay = bm >= .25 & bm <  .50 ;
        good = bm >= .50 & bm <  .75 ;
        best = bm >= .75 ;
    
    case{'tria3','dual3'}
        poor = bm <  .30 ;
        okay = bm >= .30 & bm <  .60 ;
        good = bm >= .60 & bm <  .90 ;
        best = bm >= .90 ;
        
    end
 
    r = [.85,.00,.00] ; y = [1.0,.95,.00] ;
    g = [.00,.90,.00] ; k = [.60,.60,.60] ;
    
    bar(bm(poor),hc(poor),1.05,...
        'facecolor',r,'edgecolor',r) ;
    bar(bm(okay),hc(okay),1.05,...
        'facecolor',y,'edgecolor',y) ;
    bar(bm(good),hc(good),1.05,...
        'facecolor',g,'edgecolor',g) ;
    bar(bm(best),hc(best),1.05,...
        'facecolor',k,'edgecolor',k) ;
    
    axis tight;    
    set(gca,'ycolor', get(gca,'color'),'ytick',[],...
        'xtick',.0:.2:1.,'layer','top','fontsize',...
            14,'linewidth',2.,'ticklength',[.025,.025],...
                'box','off','xlim',[0.,1.]) ;
    
    mins = max(0.010,min(sc)); %%!! so that axes don't obscure!
    maxs = min(0.990,max(sc));
    
    line([ mins, mins],...
        [0,max(hc)],'color','r','linewidth',1.5);
    line([mean(sc),mean(sc)],...
        [0,max(hc)],'color','r','linewidth',1.5);
    
    if ( mins > .4)
        text(mins-.01,.9*max(hc),num2str(min(sc),'%16.3f'),...
            'horizontalalignment',...
                'right','fontsize',15) ;
    else
        text(mins+.01,.9*max(hc),num2str(min(sc),'%16.3f'),...
            'horizontalalignment',...
                'left' ,'fontsize',15) ;
    end
    
    if ( mean(sc) > mins + .150)
    text(mean(sc)-.01,.9*max(hc),num2str(mean(sc),'%16.3f'),...
        'horizontalalignment','right','fontsize',15) ;
    end
    
    switch (ty)
    case 'tria4'
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-.04,0.0, ...
        '$\mathcal{Q}^{\mathcal{T}}_{\tau}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-.04,0.0,'Q^{t}_{\tau}',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
        
    case 'tria3'
        
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-.04,0.0, ...
        '$\mathcal{Q}^{\mathcal{T}}_{f}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-.04,0.0,'Q^{t}_{f}',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
           
    case 'dual4'
        
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-.04,0.0, ...
        '$\mathcal{Q}^{\mathcal{D}}_{\tau}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-.04,0.0,'Q^{d}_{\tau}',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
        
    case 'dual3'
        
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
        text(-.04,0.0, ...
        '$\mathcal{Q}^{\mathcal{D}}_{f}$',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter','latex') ;
    else
        text(-.04,0.0,'Q^{d}_{f}',...
            'horizontalalignment','right',...
                'fontsize',22,'interpreter',  'tex') ;
    end
        
    end
    
end

function hfnhist(hf,ty)
%HFNHIST draw histogram for "hfunc" quality-metric.

    be = linspace(0.,2.,101);
    bm = (be(1:end-1)+be(2:end)) / 2.;
    hc = histc(hf,be);

    poor = bm <  .40 | bm >= 1.6  ;
    okay =(bm >= .40 & bm <  .60 )| ...
          (bm >= 1.4 & bm <  1.6 );
    good =(bm >= .60 & bm <  .80 )| ...
          (bm >= 1.2 & bm <  1.4 );
    best = bm >= .80 & bm <  1.2 ;
 
    r = [.85,.00,.00] ; y = [1.0,.95,.00] ;
    g = [.00,.90,.00] ; k = [.60,.60,.60] ;
    
    bar(bm(poor),hc(poor),1.05,...
        'facecolor',r,'edgecolor',r) ;
    bar(bm(okay),hc(okay),1.05,...
        'facecolor',y,'edgecolor',y) ;
    bar(bm(good),hc(good),1.05,...
        'facecolor',g,'edgecolor',g) ;
    bar(bm(best),hc(best),1.05,...
        'facecolor',k,'edgecolor',k) ;
    
    axis tight; 
    set(gca,'ycolor', get(gca,'color'),'ytick',[],...
        'xtick',.0:.5:2.,'layer','top','fontsize',...
            14,'linewidth',2.,'ticklength',[.025,.025],...
                'box','off','xlim',[0.,2.]);
    
    line([max(hf),max(hf)],...
        [0,max(hc)],'color','r','linewidth',1.5);
    
    text(max(hf)+.02,.90*max(hc),num2str(max(hf),'%16.2f'),...
        'horizontalalignment','left','fontsize',15) ;
        
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
    
    text(max(hf)-.18,.45*max(hc),'$\bar{\sigma}_{h}\! = $',...
        'horizontalalignment','left',...
            'fontsize',16,'interpreter','latex') ;

    text(max(hf)+.02,.45*max(hc),num2str(mad(hf),'%16.2f'),...
        'horizontalalignment','left','fontsize',15) ;
        
    end
    
    if ( ~(exist('OCTAVE_VERSION','builtin') > +0) )
    
    text(-0.10,0.0,'$h_{r}$','horizontalalignment','right',...
        'fontsize',22,'interpreter','latex') ;

    else
    
    text(-0.10,0.0, 'h_{r}' ,'horizontalalignment','right',...
        'fontsize',22,'interpreter',  'tex') ;
    
    end
    
end



