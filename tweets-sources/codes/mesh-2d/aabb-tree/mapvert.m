function [tm,im] = mapvert(tr,pi)
%MAPVERT find the tree-to-vertex mappings.
%   [TM,IM] = MAPVERT(TR,PI) returns the tree-to-vertex and 
%   vertex-to-tree mappings for a given aabb-tree TR and a 
%   collection of query vertices PI.
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
%   See also QUERYSET, MAPRECT, MAKETREE

%   Darren Engwirda : 2014 --
%   Email           : de2363@columbia.edu
%   Last updated    : 06/04/2017

%----------------------- call SCANTREE to do the actual work
    if (nargout == +1)
       [tm   ] = scantree(tr,pi,@partvert);     
    else
       [tm,im] = scantree(tr,pi,@partvert);
    end
    
end

function [j1,j2] = partvert(pi,b1,b2)
%PARTVERT partition points between boxes B1,B2 for SCANTREE.

    j1 = true(size(pi,1),1);
    j2 = true(size(pi,1),1);

    nd = size(b1,2) / +2;
    
    for ax = +1 : nd
%--------------- remains TRUE if inside bounds along axis AX
    j1 = j1 & pi(:,ax)>=b1(ax+nd*0) ...
            & pi(:,ax)<=b1(ax+nd*1) ;
        
%--------------- remains TRUE if inside bounds along axis AX
    j2 = j2 & pi(:,ax)>=b2(ax+nd*0) ...
            & pi(:,ax)<=b2(ax+nd*1) ;    
    end

end



