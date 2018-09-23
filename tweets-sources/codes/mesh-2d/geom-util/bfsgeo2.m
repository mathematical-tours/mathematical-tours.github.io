function [node,PSLG,part] = bfsgeo2(node,PSLG,seed)
%BFSGEO2 partition geometry about "seeds" via breadth-first
%search.
%   [NODE,EDGE,PART] = BFSGEO2(NODE,EDGE,SEED) returns a set
%   of 2-manifold geometry partitions by expanding about a
%   list of "seed" points. Partitions expand in a breadth-
%   first sense until a polygon bounday is encountered. SEED
%   is an S-by-2 array of XY-coordinates to expand around,
%   NODE is an N-by-2 array of polygon vertices, and EDGE is
%   an E-by-2 array of polygon edge indexing. Each row in 
%   EDGE represents an edge of the polygon, such that
%   NODE(EDGE(JJ,1),:) and NODE(EDGE(JJ,2),:) are the XY-co-
%   ordinates of the endpoints of the JJ-TH edge. PART is an
%   S-by-1 cell array of geometry partitions, where each
%   PART{KK} is a list of edge indices into EDGE that define 
%   the KK-TH partition.
%
%   This function may be useful when seeking to partition
%   complex, non-manifold geometry into a format that's app-
%   ropriate for the triangulation algorithms in REFINE2. 
%
%   See also REFINE2, FIXGEO2, BFSTRI2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 11/10/2017
%-----------------------------------------------------------

%---------------------------------------------- basic checks    
    if ( ~isnumeric(node) || ...
         ~isnumeric(PSLG) || ...
         ~isnumeric(seed) )
        error('bfsgeo2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(node) ~= +2 || ...
        ndims(PSLG) ~= +2 || ...
        ndims(seed) ~= +2 )
        error('bfsgeo2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (size(node,2)~= +2 || ...
        size(PSLG,2)~= +2 || ...
        size(seed,2)~= +2 )
        error('bfsgeo2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    nnod = size(node,1) ;
    nedg = size(PSLG,1) ;
    
%---------------------------------------------- basic checks
    if (min([PSLG(:)])<+1 || ...
        max([PSLG(:)])>nnod)
        error('bfsgeo2:invalidInputs', ...
            'Invalid EDGE input array.') ;
    end
    
%------------------------------------------ assemble full CDT
   [node,PSLG,tria] = deltri2 (node,PSLG) ;

%------------------------------------------ find seeds in CDT
   [sptr,stri] = findtria(node,tria,seed) ;

    okay = sptr(:,2) ...
        >= sptr(:,1) ;
    itri = stri(sptr(okay,1));

%------------------------------------------ PART for all seed
    part = {} ;

    for ipos = +1 : size (itri,1)

    %-------------- BFS about current tria.
       [mark] = ...
            bfstri2(PSLG,tria,itri(ipos)) ;
   
    %-------------- match tria./poly. edges
        edge = [
            tria(mark,[1,2]) ;
            tria(mark,[2,3]) ;
            tria(mark,[3,1]) ;
               ] ;
           
        edge = sort(edge,+2) ;
        PSLG = sort(PSLG,+2) ;
      
       [same,epos] = ...
            setset2(edge(:,1:2),PSLG) ;
     
    %-------------- find match multiplicity
        epos = epos(epos>+0) ;
        epos = sort(epos);
       
        eidx = ...
            find(diff(epos)) ;
        
        eptr = ...
            [1;eidx+1;length(epos)+1] ;
        enum = ...
            eptr(2:end)-eptr(1:end-1) ;
 
    %---------- select singly-matched edges
        part{ipos} = ...
            epos(eptr(enum == 1)) ;
       
    end

end



