function [mesh] = loadmsh(name)
%LOADMSH load a *.MSH file for JIGSAW.
%
%   MESH = LOADMSH(NAME);
%
%   The following are optionally read from "NAME.MSH". Enti-
%   ties are loaded if they are present in the file:
%
%   .IF. MESH.MSHID == 'EUCLIDEAN-MESH':
%   -----------------------------------
%
%   MESH.POINT.COORD - [NPxND+1] array of point coordinates, 
%       where ND is the number of spatial dimenions. 
%       COORD(K,ND+1) is an ID tag for the K-TH point.
%
%   MESH.POINT.POWER - [NPx 1] array of vertex "weights", 
%       associated with the dual "power" tessellation.
%
%   MESH.EDGE2.INDEX - [N2x 3] array of indexing for EDGE-2 
%       elements, where INDEX(K,1:2) is an array of 
%       "point-indices" associated with the K-TH edge, and 
%       INDEX(K,3) is an ID tag for the K-TH edge.
%
%   MESH.TRIA3.INDEX - [N3x 4] array of indexing for TRIA-3 
%       elements, where INDEX(K,1:3) is an array of 
%       "point-indices" associated with the K-TH tria, and 
%       INDEX(K,4) is an ID tag for the K-TH tria.
%
%   MESH.QUAD4.INDEX - [N4x 5] array of indexing for QUAD-4 
%       elements, where INDEX(K,1:4) is an array of 
%       "point-indices" associated with the K-TH quad, and 
%       INDEX(K,5) is an ID tag for the K-TH quad.
%
%   MESH.TRIA4.INDEX - [M4x 5] array of indexing for TRIA-4 
%       elements, where INDEX(K,1:4) is an array of 
%       "point-indices" associated with the K-TH tria, and 
%       INDEX(K,5) is an ID tag for the K-TH tria.
%
%   MESH.HEXA8.INDEX - [M8x 9] array of indexing for HEXA-8 
%       elements, where INDEX(K,1:8) is an array of 
%       "point-indices" associated with the K-TH hexa, and 
%       INDEX(K,9) is an ID tag for the K-TH hexa.
%
%   MESH.WEDG6.INDEX - [M6x 7] array of indexing for WEDG-6 
%       elements, where INDEX(K,1:6) is an array of 
%       "point-indices" associated with the K-TH wedg, and 
%       INDEX(K,7) is an ID tag for the K-TH wedg.
%
%   MESH.PYRA5.INDEX - [M5x 6] array of indexing for PYRA-5 
%       elements, where INDEX(K,1:5) is an array of 
%       "point-indices" associated with the K-TH pyra, and 
%       INDEX(K,6) is an ID tag for the K-TH pyra.
%
%   MESH.VALUE - [NPxNV] array of "values" associated with
%       the vertices of the mesh.
%
%
%   .IF. MESH.MSHID == 'ELLIPSOID-MESH':
%   -----------------------------------
%
%   MESH.RADII - [ 3x 1] array of principle ellipsoid radii.
%
%
%   .IF. MESH.MSHID == 'EUCLIDEAN-GRID':
%   .OR. MESH.MSHID == 'ELLIPSOID-GRID':
%   -----------------------------------
%
%   MESH.POINT.COORD - [NDx1] cell array of grid coordinates 
%       where ND is the number of spatial dimenions. Each
%       array COORD{ID} should be a vector of grid coord.'s,
%       increasing or decreasing monotonically.
%
%   MESH.VALUE - [NMxNV] array of "values" associated with 
%       the vertices of the grid, where NM is the product of
%       the dimensions of the grid. NV values are associated 
%       with each vertex.
%
%   See also JIGSAW, SAVEMSH
%

%-----------------------------------------------------------
%   Darren Engwirda
%   github.com/dengwirda/jigsaw/
%   11-Nov-2017
%   darren.engwirda@columbia.edu
%-----------------------------------------------------------
%

    mesh = [] ;

    try

    ffid = fopen(name,'r') ;
    
    real = '%f;' ;
    ints = '%i;' ;
    
    kind = 'EUCLIDEAN-MESH';
    
    nver = +0 ;
    ndim = +0 ;
    
    while (true)
  
    %-- read next line from file
        lstr = fgetl(ffid);
        
        if (ischar(lstr) )
        
        if (length(lstr) > +0 && lstr(1) ~= '#')

        %-- tokenise line about '=' character
            tstr = regexp(lower(lstr),'=','split');
           
            switch (strtrim(tstr{1}))
            case 'mshid'

        %-- read "MSHID" data
        
                stag = regexp(tstr{2},';','split');
        
                nver = str2double(stag{1}) ;
                
                if (length(stag) >= +2)
                    kind = ...
                    strtrim(upper(stag{2}));
                end

            case 'ndims'

        %-- read "NDIMS" data
        
                ndim = str2double(tstr{2}) ;

            case 'point'

        %-- read "POINT" data

                nnum = str2double(tstr{2}) ;

                numr = nnum*(ndim+1);

                data = ...
            fscanf(ffid,[repmat(real,1,ndim),'%i'],numr);

                if (ndim == +2)
                mesh.point.coord = [ ...
                    data(1:3:end), ...
                    data(2:3:end), ...
                    data(3:3:end)] ;
                end
                if (ndim == +3)
                mesh.point.coord = [ ...
                    data(1:4:end), ...
                    data(2:4:end), ...
                    data(3:4:end), ...
                    data(4:4:end)] ;
                end
                
            case 'coord'

        %-- read "coord" data
              
                stag = regexp(tstr{2},';','split');
        
                idim = str2double(stag{1}) ;
                cnum = str2double(stag{2}) ;

                data = fscanf(ffid,'%f',cnum) ;
        
                mesh.point.coord{idim}= data  ;
                
            case 'edge2'

        %-- read "EDGE2" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 3;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,2),'%i'],numr);
                
                mesh.edge2.index = [ ...
                    data(1:3:end), ...
                    data(2:3:end), ...
                    data(3:3:end)] ;
            
                mesh.edge2.index(:,1:2) = ...
                mesh.edge2.index(:,1:2) + 1;
                
            case 'tria3'

        %-- read "TRIA3" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 4;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,3),'%i'],numr);
                
                mesh.tria3.index = [ ...
                    data(1:4:end), ...
                    data(2:4:end), ...
                    data(3:4:end), ...
                    data(4:4:end)] ;
            
                mesh.tria3.index(:,1:3) = ...
                mesh.tria3.index(:,1:3) + 1;
            
            case 'quad4'

        %-- read "QUAD4" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 5;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,4),'%i'],numr);
                
                mesh.quad4.index = [ ...
                    data(1:5:end), ...
                    data(2:5:end), ...
                    data(3:5:end), ...
                    data(4:5:end), ...
                    data(5:5:end)] ;
            
                mesh.quad4.index(:,1:4) = ...
                mesh.quad4.index(:,1:4) + 1;
        
            case 'tria4'

        %-- read "TRIA4" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 5;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,4),'%i'],numr);
                
                mesh.tria4.index = [ ...
                    data(1:5:end), ...
                    data(2:5:end), ...
                    data(3:5:end), ...
                    data(4:5:end), ...
                    data(5:5:end)] ;
            
                mesh.tria4.index(:,1:4) = ...
                mesh.tria4.index(:,1:4) + 1;
            
            case 'hexa8'

        %-- read "HEXA8" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 9;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,8),'%i'],numr);
                
                mesh.hexa8.index = [ ...
                    data(1:9:end), ...
                    data(2:9:end), ...
                    data(3:9:end), ...
                    data(4:9:end), ...
                    data(5:9:end), ...
                    data(6:9:end), ...
                    data(7:9:end), ...
                    data(8:9:end), ...
                    data(9:9:end)] ;
            
                mesh.hexa8.index(:,1:8) = ...
                mesh.hexa8.index(:,1:8) + 1;
            
            case 'wedg6'

        %-- read "WEDG6" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 7;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,6),'%i'],numr);
                
                mesh.wedg6.index = [ ...
                    data(1:7:end), ...
                    data(2:7:end), ...
                    data(3:7:end), ...
                    data(4:7:end), ...
                    data(5:7:end), ...
                    data(6:7:end), ...
                    data(7:7:end)] ;
            
                mesh.wedg6.index(:,1:6) = ...
                mesh.wedg6.index(:,1:6) + 1;
     
            case 'pyra5'

        %-- read "PYRA5" data

                nnum = str2double(tstr{2}) ;

                numr = nnum * 6;
                
                data = ...
            fscanf(ffid,[repmat(ints,1,5),'%i'],numr);
                
                mesh.pyra5.index = [ ...
                    data(1:6:end), ...
                    data(2:6:end), ...
                    data(3:6:end), ...
                    data(4:6:end), ...
                    data(5:6:end), ...
                    data(6:6:end)] ;
            
                mesh.pyra5.index(:,1:5) = ...
                mesh.pyra5.index(:,1:5) + 1;
            
            case 'value'

        %-- read "VALUE" data

                stag = regexp(tstr{2},';','split');
                
                nnum = str2double(stag{1}) ;
                vnum = str2double(stag{2}) ;
               
                numr = nnum * (vnum+1) ;
                
                fstr = repmat(real,1,vnum) ;
                
                data = fscanf( ...
                  ffid,fstr(1:end-1),numr) ;
                
                mesh.value = zeros(nnum, vnum);
                
                for vpos = +1 : vnum
                
                mesh.value(:,vpos) = ...
                    data(vpos:vnum:end) ;

                end
               
            case 'power'

        %-- read "POWER" data

                stag = regexp(tstr{2},';','split');
                
                nnum = str2double(stag{1}) ;
                pnum = str2double(stag{2}) ;
               
                numr = nnum * (pnum+0) ;
                
                fstr = repmat(real,1,pnum) ;
                
                data = fscanf( ...
                  ffid,fstr(1:end-1),numr) ;
                
                mesh.point.power = ...
                    zeros(nnum, pnum) ;
                
                for ppos = +1 : pnum
                
                mesh.point.power(:,ppos) = ...
                    data(ppos:pnum:end) ;

                end
                
            end
                       
        end
           
        else
    %-- if(~ischar(lstr)) //i.e. end-of-file
            break ;
        end
        
    end
    
    mesh.mshID = kind ;
    mesh.fileV = nver ;
    
    fclose(ffid) ;
    
    catch err

%-- ensure that we close the file regardless!
    if (ffid>-1)
    fclose(ffid) ;
    end
    
    rethrow(err) ;
    
    end

end


