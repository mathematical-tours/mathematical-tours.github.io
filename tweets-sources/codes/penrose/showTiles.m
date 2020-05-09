function showTiles(T,acolor,bcolor)
%showTiles Show Penrose rhombus tiles.
%
%   showTiles(T,acolor,bcolor) displays the Penrose rhombus tiles
%   constructed from the triangles in the input table, T. acolor is the
%   color of the rhombus tiles formed from A and A' triangles. bcolor is
%   the color of the rhombus tiles formed from B and B' triangles. If
%   not specified, bcolor is a shade of orange, and acolor is a shade of
%   blue.
%
%   EXAMPLE
%   Decompose a B triangle 4 times and display the resulting rhombus
%   tiles.
%
%       t = bTriangle([],-1,1);
%       for k = 1:4
%           t = decomposeTriangles(t);
%       end
%       showTiles(t)

%   Copyright 2018 The MathWorks, Inc.

ax = gca;

if nargin < 3
   bcolor = [253 174 97]/255;
end
if nargin < 2
   acolor = [171 217 233]/255;
end

% Display all the tiles as one patch object, constructed in face-vertex
% form.

points = [T.Apex T.Left T.Right];
points = points.';
points = points(:);
vertices = [real(points) imag(points)];
faces = reshape(1:(3*height(T)),3,[])';
colors = zeros(height(T),3);

for k = 1:height(T)
    switch T.Type(k)
        case {'A' 'Ap'}
            colors(k,:) = acolor;
            
        case {'B' 'Bp'}
            colors(k,:) = bcolor;
    end
end

patch(ax,'Faces',faces,'Vertices',vertices,...
    'EdgeColor','none',...
    'FaceColor','flat',...
    'FaceVertexCData',colors);
 
% Now display the rhombus edges. These are the triangle sides, omitting
% the triangle bases.

points = [T.Left T.Apex T.Right];
points(:,end+1) = NaN;
points = points.';
points = points(:);
x = real(points);
y = imag(points);
hold on
plot(ax,x,y,'LineWidth',1)
hold off

axis(ax,'equal')
axis(ax,'off')