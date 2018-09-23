function [vert,conn,tria,tnum] = tridiv2(varargin)
%TRIDIV2 "quadtree" refinement for 2-simplex triangulations.
%   [VERT,EDGE,TRIA,TNUM] = TRIDIV2(VERT,EDGE,TRIA,TNUM) re-
%   turns a globally refined triangulation, in which all ed-
%   ges are bisected about their midpoints. Such refinement
%   splits each triangle into four new sub-triangles accord-
%   ing to a shape-preserving scheme.
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
%   [VERT,EDGE,TRIA,TNUM] = TRIDIV2(... ,TDIV) returns a se-
%   lectively refined mesh, where TDIV is a T-by-1 logical 
%   array, with TDIV(KK) = TRUE if TRIA(KK,:) is to be refi-
%   ned. Such triangles are refined using the four-way split
%   described above. Additionally, a "halo" of adjacent tri-
%   angles are also refined to preseve compatibility of the
%   mesh. Such triangles are refined using a two-way bisect-
%   ion type scheme. 

%   See also REFINE2, SMOOTH2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 29/01/2017

%---------------------------------------------- extract args
    vert = []; conn = []; tria = []; tnum = [] ;
    tdiv = []; 
    
    if (nargin>=+1), vert = varargin{1}; end
    if (nargin>=+2), conn = varargin{2}; end
    if (nargin>=+3), tria = varargin{3}; end
    if (nargin>=+4), tnum = varargin{4}; end
    if (nargin>=+5), tdiv = varargin{5}; end
    
%---------------------------------------------- default arg.
    if (isempty(tnum))
        tnum = ones(size(tria,1),1) ;
    end
    if (isempty(tdiv))
        tdiv = true(size(tria,1),1) ;
    end
    
%---------------------------------------------- basic checks    
    if ( ~isnumeric(vert) || ...
         ~isnumeric(conn) || ...
         ~isnumeric(tria) || ...
         ~isnumeric(tnum) || ...
         ~islogical(tdiv) )
        error('tridiv2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    tnum = tnum(:) ; tdiv = tdiv(:) ;
    if (ndims(vert) ~= +2 || ...
        ndims(conn) ~= +2 || ...
        ndims(tria) ~= +2 || ...
        ndims(tnum) ~= +2 || ...
        ndims(tdiv) ~= +2 )
        error('tridiv2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || ...
        size(conn,2)~= +2 || ...
        size(tria,2)~= +3 || ...
        size(tnum,2)~= +1 || ...
        size(tria,1)~= size(tnum,1) )
        error('tridiv2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    nvrt = size(vert,1) ;

%---------------------------------------------- basic checks
    if (min(min(conn(:,1:2))) < +1 || ...
            max(max(conn(:,1:2))) > nvrt )
        error('tridiv2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
    if (min(min(tria(:,1:3))) < +1 || ...
            max(max(tria(:,1:3))) > nvrt )
        error('tridiv2:invalidInputs', ...
            'Invalid TRIA input array.') ;
    end

%------------------------------ assemble extended adj. info.
   [edge,tria] = tricon2(tria,conn) ;

    ediv = false(size(edge,1),1) ;
    ediv(tria(tdiv,4:6))  = true ;

    snum = length(find(ediv));

    while (true)

    %-------------------------- tria's with >= 2 edge splits
        div3 = sum(double( ...
            ediv(tria(:,4:6))),2)>=2;
    
    %-------------------------- expand onto adj. edge splits
        ediv(tria(div3,4:6)) = true ; 
    
        snew = length(find(ediv)) ;
 
        if (snew == snum), break; end
 
        snum = snew ;
    
    end

%------------------------------ tria's with == 1 edge splits
    div1 = sum( ...
    double( ediv(tria(:,4:6))),2)==1;

%------------------------------ indexing for mid-point vert.
    ivec = zeros(size(edge,1),1);
    ivec(ediv) = ...
        (1:snum)' + size(vert,1);

%------------------------------------------ update vert. set
    emid = vert(edge(ediv,1),:) ...
         + vert(edge(ediv,2),:) ;
    vert =[vert ; emid * 0.5] ;
    
%------------------------------------------ update edge. set
    [cvec,eloc] = ...
        setset2(conn,edge(ediv, 1:2)) ;
 
    epos = ivec(ediv);
    epos = epos(eloc(eloc>0)) ;
 
    conn = [conn(~cvec,1:2) ; ...
        conn(cvec,1) , epos ; ...
        conn(cvec,2) , epos ] ;
    
%------------------------------------ push 1-to-4 refinement  
    tr41 = ones(length(find(div3)),3) ;
    tn41 = ones(length(find(div3)),1) ;
    tn41(:,1) = tnum(div3,1);
    tr41(:,1) = tria(div3,1);
    tr41(:,2) = ivec(tria(div3,4));
    tr41(:,3) = ivec(tria(div3,6));
    
    tr42 = ones(length(find(div3)),3) ;
    tn42 = ones(length(find(div3)),1) ;
    tn42(:,1) = tnum(div3,1);
    tr42(:,1) = ivec(tria(div3,4));
    tr42(:,2) = tria(div3,2);
    tr42(:,3) = ivec(tria(div3,5));
    
    tr43 = ones(length(find(div3)),3) ;
    tn43 = ones(length(find(div3)),1) ;
    tn43(:,1) = tnum(div3,1);
    tr43(:,1) = ivec(tria(div3,6));
    tr43(:,2) = ivec(tria(div3,5));
    tr43(:,3) = tria(div3,3);
    
    tr44 = ones(length(find(div3)),3) ;
    tn44 = ones(length(find(div3)),1) ;
    tn44(:,1) = tnum(div3,1);
    tr44(:,1) = ivec(tria(div3,6));
    tr44(:,2) = ivec(tria(div3,4));
    tr44(:,3) = ivec(tria(div3,5));
    
%----------------------- push 1-to-2 refinement about edge 1
    tvec = false(size(tria,1), 1) ;
    tvec(ediv(tria(:,4))&div1) = true ;
    
    tr21 = ones(length(find(tvec)),3) ;
    tn21 = ones(length(find(tvec)),1) ;
    tn21(:,1) = tnum(tvec,1);
    tr21(:,1) = ivec(tria(tvec,4));
    tr21(:,2) = tria(tvec,3);
    tr21(:,3) = tria(tvec,1);
    
    tr22 = ones(length(find(tvec)),3) ;
    tn22 = ones(length(find(tvec)),1) ;
    tn22(:,1) = tnum(tvec,1);
    tr22(:,1) = ivec(tria(tvec,4));
    tr22(:,2) = tria(tvec,2);
    tr22(:,3) = tria(tvec,3);

%----------------------- push 1-to-2 refinement about edge 2
    tvec = false(size(tria,1), 1) ;
    tvec(ediv(tria(:,5))&div1) = true ;
    
    tr23 = ones(length(find(tvec)),3) ;
    tn23 = ones(length(find(tvec)),1) ;
    tn23(:,1) = tnum(tvec,1);
    tr23(:,1) = ivec(tria(tvec,5));
    tr23(:,2) = tria(tvec,1);
    tr23(:,3) = tria(tvec,2);
    
    tr24 = ones(length(find(tvec)),3) ;
    tn24 = ones(length(find(tvec)),1) ;
    tn24(:,1) = tnum(tvec,1);
    tr24(:,1) = ivec(tria(tvec,5));
    tr24(:,2) = tria(tvec,3);
    tr24(:,3) = tria(tvec,1);
    
%----------------------- push 1-to-2 refinement about edge 3
    tvec = false(size(tria,1), 1) ;
    tvec(ediv(tria(:,6))&div1) = true ;
    
    tr25 = ones(length(find(tvec)),3) ;
    tn25 = ones(length(find(tvec)),1) ;
    tn25(:,1) = tnum(tvec,1);
    tr25(:,1) = ivec(tria(tvec,6));
    tr25(:,2) = tria(tvec,2);
    tr25(:,3) = tria(tvec,3);
    
    tr26 = ones(length(find(tvec)),3) ;
    tn26 = ones(length(find(tvec)),1) ;
    tn26(:,1) = tnum(tvec,1);
    tr26(:,1) = ivec(tria(tvec,6));
    tr26(:,2) = tria(tvec,1);
    tr26(:,3) = tria(tvec,2);
    
%------------------------------------------ update tria. set
    tria = [tria(~div1&~div3,1:3) ; ...
        tr41 ; tr42 ; tr43 ; tr44 ; ...
        tr21 ; tr22 ; 
        tr23 ; tr24 ; 
        tr25 ; tr26 ] ;
        
    tnum = [tnum(~div1&~div3,1:1) ; ...
        tn41 ; tn42 ; tn43 ; tn44 ; ...
        tn21 ; tn22 ; 
        tn23 ; tn24 ; 
        tn25 ; tn26 ] ;
   
end



