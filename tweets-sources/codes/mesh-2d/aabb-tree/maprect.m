function [tm,im] = maprect(tr,pr)
%MAPRECT find the tree-to-rectangle mappings.
%   [TM,IM] = MAPRECT(TR,PR) returns the tree-to-rectangle 
%   and rectangle-to-tree mappings for a given aabb-tree TR 
%   and a collection of query vertices PI.
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
%   See also QUERYSET, MAPVERT, MAKETREE

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 09/04/2017

%----------------------- call SCANTREE to do the actual work
    if (nargout == +1)
       [tm   ] = scantree(tr,pr,@partrect);     
    else
       [tm,im] = scantree(tr,pr,@partrect);
    end
    
end

function [j1,j2] = partrect(pr,b1,b2)
%PARTRECT partition points between boxes B1,B2 for SCANTREE.

    j1 = true(size(pr,1),1) ;
    j2 = true(size(pr,1),1) ;

    nd = size(b1,2) / +2;
    
    for ax = +1 : nd
%--------------- remains TRUE if inside bounds along axis AX
    j1 = j1 & pr(:,ax+nd*1) ...
           >= b1(  ax+nd*0) ...
            & pr(:,ax+nd*0) ...
           <= b1(  ax+nd*1) ;
        
%--------------- remains TRUE if inside bounds along axis AX
    j2 = j2 & pr(:,ax+nd*1) ...
           >= b2(  ax+nd*0) ...
            & pr(:,ax+nd*0) ...
           <= b2(  ax+nd*1) ;    
    end

end


