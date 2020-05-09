function Vcells = clippedVoronoi(nodes,area)
%
% function Vcells = clippedvoronoi(nodes, area)
% nodes  = [Xnodes Ynodes]
% area = [Xarea Yarea]
% Clipped Voronoi Tesselation - Version 1.0
%
%   clippedVoronoi(nodes,area) plots the  Voronoi diagram
%   derived from the nodes with coordinates nodes = [Xnodes Ynodes], while bounding
%   the result in the area described by the vertices area = [Xarea Yarea]. Vcells
%   are the vertices describing the voronoi cells. 
%
%
% NOTE: This function makes use of the voronoi() function of Matlab as well 
% as the polybool() function to intersect the diagram with the area.
%
% Example:
%   nodes =[1 2; 2 1.5; 2.5 2.5;1.5 2.7;2.7 1.7];
%   area =[ 1 3; 3 3; 3 1; 1 1; 0.5 2];
%   Vcells = clippedVoronoi(nodes,area)
%
%
%
%                                                               
%  Yiannis Kantaros, Alexandros Osana, July 2016
%


N = size(nodes,1);% Number of sensors

clear nodesExtended
nodesExtended = nodes;

nodesExtended(N+1,2)=-10000;...needed for Voronoi computation of Infinite Cells
nodesExtended(N+2,2)=-10000;...%...
nodesExtended(N+3,2)=10000;...%...
nodesExtended(N+4,2)=10000;...
nodesExtended(N+1,1)=-10000;...
nodesExtended(N+2,1)=10000;...
nodesExtended(N+3,1)=-10000;...
nodesExtended(N+4,1)=10000;...


[vertices,cells] = voronoin([nodesExtended(:,1),nodesExtended(:,2)]); %Performs tesselation and returns vertices and cells (as Indices of verticed) of the Tesselation

% ========== Plotting ============
figure()
hold on

%  Voronoi Vertices 
inAreaFlag = inpolygon(vertices(:,1),vertices(:,2),area(:,1),area(:,2));%Plotting vertices only in the area of the vineyard
for i = 1:1:length(vertices)
    if inAreaFlag(i) == 1
        plot(vertices(i,1),vertices(i,2),'y+','LineWidth',2)    
    end
end


% Area and nodes
for k = 1:size(area,1)-1
   plot(area(k:k+1,1),area(k:k+1,2),'k-','LineWidth',3) 
end
plot([area(end,1) area(1,1)], [area(end,2) area(1,2)],'k-','LineWidth',3)
plot(nodes(1:N,1),nodes(1:N,2),'r.')
for i = 1:N; text(nodes(i,1),nodes(i,2),['C' num2str(i) ' '],'HorizontalAlignment','right');end % Ploting number of node (same as cell)

%Cells 
for e=1:1:N
    vorCellX=vertices(cells{e},1);
    vorCellY=vertices(cells{e},2);
    
    [clippedCellX,clippedCellY]=polybool('&',area(:,1),area(:,2),vorCellX,vorCellY);%Intersection of cell with the area (Clipped Cell)
    
    Vcells{e,1}=clippedCellX;
    Vcells{e,2}=clippedCellY;
    
    plot(cell2mat(Vcells(e,1)),cell2mat(Vcells(e,2)),'b-')
end
axis equal
hold off;

return 