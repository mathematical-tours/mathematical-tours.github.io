function [ffun,flag] = limgrad(edge,elen,ffun,dfdx,imax)
%LIMGRAD impose "gradient-limits" on a function defined over 
%an undirected graph.
%   [FNEW] = LIMGRAD(EDGE,ELEN,FFUN,DFDX,ITER) computes a
%   "gradient-limited" function FNEW on the undirected graph
%   {EDGE,ELEN}, where EDGE is an NE-by-2 array of edge ind-
%   ices, and ELEN is an NE-by-1 array of edge lengths. 
%   Gradients are limited over the graph edges, such that
%
%       ABS(FNEW(N2)-FNEW(N1)) <= ELEN(II) * DFDX,
%
%   where N1=EDGE(II,1) and N2=EDGE(II,2) are the two nodes 
%   in the II-TH edge. An iterative algorithm is used, swee-
%   ping over an "active-set" of graph edges until converge-
%   nce is achieved. A maximum of IMAX iterations are done.
%
%   [FNEW,FLAG] = LIMGRAD(...) also returns a boolean FLAG,
%   with FLAG=TRUE denoting convergence. 
%
%   See also LIMHFN2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 18/04/2017

%---------------------------------------------- basic checks    
    if ( ~isnumeric(edge) || ...
         ~isnumeric(elen) || ...
         ~isnumeric(ffun) || ...
         ~isnumeric(dfdx) || ...
         ~isnumeric(imax) )
        error('limgrad:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
%---------------------------------------------- basic checks
    if (ndims(edge) ~= +2 || ...
        ndims(elen)  > +2 || ...
        ndims(ffun)  > +2 || ...
        numel(dfdx) ~= +1 || ...
        numel(imax) ~= +1 )
        error('limgrad:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(edge,2) < +2 || ...
        size(elen,2)~= +1 || ...
        size(ffun,2)~= +1 || ...
        size(edge,1)~= size(elen,1) )
        error('limgrad:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    nnod = size(ffun,1) ;

%---------------------------------------------- basic checks
    if (dfdx < +0. || imax < +0)
        error('limgrad:invalidInputArgument', ...
            'Invalid input parameter.');
    end
    if (min(min(edge(:,1:2))) < +1 || ...
            max(max(edge(:,1:2))) > nnod )
        error('limgrad:invalidInputArgument', ...
            'Invalid EDGE input array.') ;
    end

%-- IVEC(NPTR(II,1):NPTR(II,2)) are edges adj. to II-TH node
    nvec = [edge(:,1); edge(:,2)];
    ivec = [(1:size(edge,1))'; ...
            (1:size(edge,1))'] ;

   [nvec,pidx] = sort (nvec) ;
    ivec       = ivec (pidx) ;
    
    mark = false(nnod,1) ;
    mark(edge(:,1)) = true ;
    mark(edge(:,2)) = true ;
    
    idxx = find(diff(nvec) > +0) ;
    
    nptr = zeros(nnod,2) ;
    nptr(:,2) = -1 ;
    nptr(mark,1) = [+1; idxx+1];
    nptr(mark,2) = [idxx; nnod];
    
%----------------------------- ASET=ITER if node is "active"
    aset = zeros(size(ffun,1),1) ;
    
%----------------------------- exhaustive 'til all satisfied 
    ftol = min(ffun) * sqrt(eps) ;
    
    for iter = +1 : imax
    
    %------------------------- find "active" nodes this pass
        aidx = find(aset == iter - 1) ;
        
        if (isempty(aidx)), break; end
      
    %------------------------- reorder => better convergence
       [aval,idxx] = sort(ffun(aidx)) ;
        
        aidx = aidx(idxx);
       
    %------------------------- visit adj. edges and set DFDX
        for ipos = 1 : length(aidx)
            npos = aidx(ipos) ;
            for jpos = nptr(npos,1) ...
                     : nptr(npos,2)
                
                epos = ivec(jpos,1) ;
                
                nod1 = edge(epos,1) ;
                nod2 = edge(epos,2) ;

            %----------------- calc. limits about min.-value
                if (ffun(nod1) > ffun(nod2))
                
                fun1 = ffun(nod2) ...
                     + elen(epos) * dfdx ;
                
                if (ffun(nod1) > fun1+ftol)
                    ffun(nod1) = fun1;
                    aset(nod1) = iter;
                end

                else
                
                fun2 = ffun(nod1) ...
                     + elen(epos) * dfdx ;
                    
                if (ffun(nod2) > fun2+ftol)
                    ffun(nod2) = fun2;
                    aset(nod2) = iter;
                end
                
                end
                 
            end
        end
        
    end
     
    flag = (iter < imax) ;
    
end



