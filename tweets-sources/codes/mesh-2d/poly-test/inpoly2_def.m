function [stat] = inpoly2_def(vert,node,edge)
%INPOLY2_DEF the local m-code version of the crossing-number
%test. Loop over edges; do a binary-search for the first ve-
%rtex that intersects with the edge y-range; do crossing-nu-
%mber comparisons; break when the local y-range is exceeded.

    nvrt = size (vert,1) ;
    nnod = size (node,1) ;
    nedg = size (edge,1) ;

    stat = false(nvrt,1) ;

%----------------------------------- loop over polygon edges
    for epos = +1 : size(edge,1)
    
        inod = edge(epos,1) ;
        jnod = edge(epos,2) ;

    %------------------------------- calc. edge bounding-box
        yone = node(inod,2) ;
        ytwo = node(jnod,2) ;
        xone = node(inod,1) ;
        xtwo = node(jnod,1) ;
        
        ydel = ytwo - yone;
        xdel = xtwo - xone;

        xmin = min(xone,xtwo) ;

    %------------------------------- find VERT(IPOS,2)<=YONE
        ilow = +1 ; iupp = nvrt ;
        
        while (ilow < iupp - 1)        
            imid = ilow ...
            + floor((iupp-ilow) / 2);
            
            if (vert(imid,2) < yone)
                ilow = imid ;
            else
                iupp = imid ;
            end           
        end

    %------------------------------- calc. edge-intersection
        for jpos = ilow+1 : nvrt
        
            ypos = vert(jpos,2) ;
            if (ypos <  ytwo)
                xpos = vert(jpos,1) ;
                if (xpos >= xmin)
                    if ( ...
                    ydel* (xpos-xone) < ...
                    xdel* (ypos-yone) )
                    
                    stat(jpos) = ...
                        ~stat(jpos) ;
                    end
                else
                    stat(jpos) = ...
                        ~stat(jpos) ;
                end
            else
                break ;            % done -- due to the sort
            end
        
        end

    end

end



