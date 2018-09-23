function [vert,edge,tria,tnum] = triread(name)
%TRIREAD read two-dimensional triangulation data from file.
%   [VERT,EDGE,TRIA,TNUM] = TRIREAD(NAME) returns the 2-sim-
%   plex triangulation {VERT,TRIA} contained in the mesh-fi-
%   le NAME.
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
%   Data is returned non-empty if it is present in the file.
%
%   See also TRISAVE 

%   This routine borrows functionality from the JIGSAW pack-
%   age: github.com/dengwirda/jigsaw-matlab.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 01/10/2017

    vert = []; edge = []; tria = []; tnum = [];

%---------------------------------------------- basic checks
    if (~ischar(name))
        error('triread:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

%----------------------------------- borrow JIGSAW I/O func!
   [mesh] = loadmsh(name) ;

%----------------------------------- extract data if present
    if (isfield(mesh,'point') && ...
        isfield(mesh.point,'coord'))
    vert = mesh.point.coord(:,1:2) ;    
    end
    if (isfield(mesh,'edge2') && ...
        isfield(mesh.edge2,'index'))
    edge = mesh.edge2.index(:,1:2) ;
    end
    if (isfield(mesh,'tria3') && ...
        isfield(mesh.tria3,'index'))
    tria = mesh.tria3.index(:,1:3) ;
    tnum = mesh.tria3.index(:,  4) ;
    end

end



