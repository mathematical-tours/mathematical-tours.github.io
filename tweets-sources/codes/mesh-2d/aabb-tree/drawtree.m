function drawtree(tr,varargin)
%DRAWTREE draw an aabb-tree generated using MAKETREE.
%   DRAWTREE(TR) draws the tree TR for cases in R^2 and R^3.
%
%   See also MAKETREE

%   Darren Engwirda : 2014 --
%   Email           : darren.engwirda@columbia.edu
%   Last updated    : 18/12/2014

%---------------------------------------------- basic checks
    if (~isa(tr,'struct') || ...
        ~isfield(tr,'xx') || ...
        ~isfield(tr,'ii') || ...
        ~isfield(tr,'ll') )
        error('drawtree:incorrectInputClass', ...
            'Incorrect aabb-tree.') ;
    end

%----------------------------------------- find "leaf" nodes
    lf = ~cellfun('isempty', tr.ll) ;
    
    fc = [.95,.95,.55] ;
    ec = [.15,.15,.15] ;
    
%-------------------------- draw all "leaf" nodes as patches
    switch (size(tr.xx,2))
        case 4
%----------------------------------------------- tree in R^2
        np = numel(find(lf));
    %------------------------------------------------- nodes
        pp = [tr.xx(lf,1),tr.xx(lf,2)
              tr.xx(lf,3),tr.xx(lf,2)
              tr.xx(lf,3),tr.xx(lf,4)
              tr.xx(lf,1),tr.xx(lf,4)
             ] ;
    %------------------------------------------------- faces
        bb = [(1:np)'+np*0,(1:np)'+np*1,...
              (1:np)'+np*2,(1:np)'+np*3
             ] ;
            
        case 6
%----------------------------------------------- tree in R^3
        np = numel(find(lf));
    %------------------------------------------------- nodes
        pp = [tr.xx(lf,1),tr.xx(lf,2),tr.xx(lf,3)
              tr.xx(lf,4),tr.xx(lf,2),tr.xx(lf,3)
              tr.xx(lf,4),tr.xx(lf,5),tr.xx(lf,3)
              tr.xx(lf,1),tr.xx(lf,5),tr.xx(lf,3)
              tr.xx(lf,1),tr.xx(lf,2),tr.xx(lf,6)
              tr.xx(lf,4),tr.xx(lf,2),tr.xx(lf,6)
              tr.xx(lf,4),tr.xx(lf,5),tr.xx(lf,6)
              tr.xx(lf,1),tr.xx(lf,5),tr.xx(lf,6)
             ] ;
    %------------------------------------------------- faces
        bb = [(1:np)'+np*0,(1:np)'+np*1,...
              (1:np)'+np*2,(1:np)'+np*3
              (1:np)'+np*4,(1:np)'+np*5,...
              (1:np)'+np*6,(1:np)'+np*7
              (1:np)'+np*0,(1:np)'+np*3,...
              (1:np)'+np*7,(1:np)'+np*4
              (1:np)'+np*3,(1:np)'+np*2,...
              (1:np)'+np*6,(1:np)'+np*7
              (1:np)'+np*2,(1:np)'+np*1,...
              (1:np)'+np*5,(1:np)'+np*6
              (1:np)'+np*1,(1:np)'+np*0,...
              (1:np)'+np*4,(1:np)'+np*5
             ] ;
            
        otherwise
%--------------------------- what to do with a tree in R^d!?
        error('scantree:unsupportedDimension', ...
            'Unsupported tree dimensionality.') ;
    end
    
    patch('faces',bb,'vertices',pp,'facecolor',fc,...
        'edgecolor',ec,'facealpha',+.2);
    
end



