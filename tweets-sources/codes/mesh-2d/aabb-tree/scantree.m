function [tm,im] = scantree(tr,pi,fn)
%SCANTREE find the tree-to-item mappings.
%   [TM,IM] = SCANTREE(TR,PI,FN) is a low-level routine that
%   returns the tree-to-item and item-to-tree mappings for a
%   given aabb-tree TR and a query collection PI. The funct-
%   ion [KI,KJ] = FN(PJ,NI,NJ) is called internally to part-
%   ition the sub-collection PJ between the aabb-tree nodes
%   NI,NJ, where:
%
%   * KI(II) = TRUE if the II-th item intersects NI, and
%   * KJ(II) = TRUE if the II-th item intersects NJ.
%
%   The tree-to-item mapping TM is a structure representing
%   the intersection of the items PI with the tree TR. TM.II 
%   is an M-by-1 array of tree indices and TM.LL is an 
%   M-by-1 cell array of item lists. Specifically, items in 
%   the list TM.LL{JJ} intersect with the node TM.II(JJ).
%
%   The item-to-tree mapping IM is a structure representing
%   the inverse mapping. IM.II is an N-by-1 array of item
%   indices and IM.LL is an N-by-1 cell array of node lists.
%   Specifically, nodes in the list IM.LL{JJ} intersect with
%   the item IM.II(JJ).
%
%   See also QUERYSET, MAPVERT, MAPRECT, MAKETREE

%   Darren Engwirda : 2014 --
%   Email           : de2363@columbia.edu
%   Last updated    : 06/04/2017

    tm.ii = [] ; tm.ll = {} ;
    im.ii = [] ; im.ll = {} ;
    
%------------------------------ quick return on empty inputs
    if (isempty(pi)), return; end
    if (isempty(tr)), return; end  
     
%---------------------------------------------- basic checks    
    if (~isa(tr, 'struct') || ...
        ~isa(pi,'numeric') || ...
        ~isa(fn,'function_handle') )
        error('scantree:incorrectInputClass', ...
            'Invalid input class.') ;
    end
%---------------------------------------------- basic checks
    if ( ~isfield(tr,'xx') || ...
         ~isfield(tr,'ii') || ...
         ~isfield(tr,'ll') )
        error('scantree:incorrectAABBstruct', ...
            'Incorrect aabb-tree.') ;
    end
    
%----------------------------------- alloc. output/workspace
    tm.ii = zeros(size(tr.ii,1),1);
    tm.ll = cell (size(tr.ii,1),1); 
    
    ss    = zeros(size(tr.ii,1),1);
    sl    = cell (size(tr.ii,1),1);
    sl{1} = (1:size(pi,1))';
    
    tf = ~cellfun('isempty',tr.ll);
    
%---------- descend tree from root, push items amongst nodes

    ss(1) = +1; ns = +1; no = +1;
    while (ns ~= +0)
    %---------------------------------- _pop node from stack
        ni = ss(ns); ns = ns - 1;
        
        if (tf(ni))
        %-- push onto tree-item mapping -- non-empty node NI 
        %-- contains items LL
            tm.ii(no) = ni ; 
            tm.ll{no} ...
                = sl{ns+1} ;
            no = no + 1 ;
        end
        
        if (tr.ii(ni,+2)~=+0)
        %--------------------- partition amongst child nodes
            c1 = ...
            tr.ii(ni,2) + 0 ;
            c2 = ...
            tr.ii(ni,2) + 1 ;
            
        %--------------------- user-defined partitions of LL
           [j1,j2] = feval( ...
                fn,pi(sl{ns+1},:), ...
                    tr.xx(c1,:),tr.xx(c2,:)) ;
                           
        %--------------------- lists of items per child node
            l1 = sl{ns+1}(j1) ; 
            l2 = sl{ns+1}(j2) ;

            if (~isempty(l1))
        %--------------------- push nonempty node onto stack
                ns = ns + 1 ; 
                ss(ns) = c1 ; 
                sl{ns} = l1 ; 
            end
            if (~isempty(l2)) 
        %--------------------- push nonempty node onto stack
                ns = ns + 1 ; 
                ss(ns) = c2 ; 
                sl{ns} = l2 ; 
            end
            
        end
        
    end
%----------------------------------------------- trim alloc.    
    tm.ii(no:end) = [];
    tm.ll(no:end) = [];

%----------------------- compute inverse map only if desired
    if (nargout==+1), return; end
    
%----------------------- accumulate pair'd tree-item matches
    ic = cell(no-1,+1); jc = tm.ll;
    for ip = +1 : no-1
        ni = tm.ii(ip);
        ic{ip} = ...
          ni * ones(length(jc{ip}),1);
    end
    ii = vertcat(ic{:}); 
    ni = size(ii,1) ;
    jj = vertcat(jc{:}); 
    nj = size(jj,1) ;
    
    if (isempty(jj)), return; end
    
    im.ll = cell (size(pi,1),1) ;
    
%---------------------------------- invert ordering via sort
   [jj,jx] =  sort(jj) ; 
    ii = ii(jx);
    jx = find(diff(jj)~=+0);
    im.ii = [jj(jx);jj(nj)];
    jx = [+0;jx;ni];
    
%----------------------- distribute single item-tree matches
    for ip = +1 : size(im.ii,1)
        im.ll{ip} = ...
            ii(jx(ip+0)+1:jx(ip+1)+0);
    end
    
end



