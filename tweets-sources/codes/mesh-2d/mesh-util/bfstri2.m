function [seen] = bfstri2(PSLG,tria,seed)
%BFSTRI2 expand about a single seed triangle via BFS. The se-
%arch terminates when constraining edges are encountered.
%SEEN(II) is TRUE if the II-TH triangle is found in the curr-
%ent expansion.
%
%   See also BFSGEO2, REFINE2, FIXGEO2

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 01/10/2017
%-----------------------------------------------------------
  
    seen = [];

%---------------------------------------------- basic checks    
    if ( ~isnumeric(tria) || ...
         ~isnumeric(seed) )
        error('bfstri2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(tria) ~= +2 )
        error('bfstri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(tria,2)~= +3 )
        error('bfstri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%---------------------------------------------- extra checks   
    if ( ~isempty  (PSLG) )
    if ( ~isnumeric(PSLG) )
        error('bfstri2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
    if (ndims(PSLG) ~= +2 )
        error('bfstri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end  
    if (size(PSLG,2)~= +2 )
        error('bfstri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    end
    
%----------------------------------------- form adj. indices
    ntri = size (tria,1);

   [edge,tria] = tricon2 (tria,PSLG);

    seed = seed(:) ;

    list = zeros(ntri,1);
    nlst = length (seed);
    list(1:nlst) = seed ;

%----------------------------------------- do BFS iterations
    seen = false(ntri,1);

    while (nlst >= +1)
        
    %-------------- pop tria from stack top
        next = list(nlst);
        nlst = nlst-1 ;    
        seen(next) = true;
    
    %-------------- visit 1-ring neighbours
        for eadj = +1 : +3
        
            epos = tria(next,eadj+3);

        %---------- find adjacent triangles
            if (edge(epos,5) == 0)
            
            if (next ~= edge(epos,3))
                tadj  = edge(epos,3);
            else
                tadj  = edge(epos,4);
            end         

            if (tadj > +0 && ~seen(tadj))
                
        %---------- add unvisited neighbour
                seen(tadj) = true ;
                nlst = nlst+1 ;
                list(nlst) = tadj ;
            
            end
                        
            end
                 
        end
        
    end

end



