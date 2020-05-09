function showDecoratedTiles(T)
%showDecoratedTiles Show Penrose rhombus tiles with connecting arcs.
%
%   showDecoratedTiles(T) displays the Penrose rhombus tiles constructed
%   from the triangles in the input table, T. Each triangle is decorated
%   with arcs so that the arcs connect smoothly from triangle to
%   triangle, resulting in an interesting geometric pattern overlaid on
%   the Penrose tiles.
%
%   EXAMPLE
%   Decompose a B triangle 4 times and display the resulting rhombus
%   tiles.
%
%       t = bTriangle([],-1,1);
%       for k = 1:4
%           t = decomposeTriangles(t);
%       end
%       showDecoratedTiles(t)

%   Copyright 2018 The MathWorks, Inc.

showTiles(T);

[arc1,arc2] = triangleCurves(T);

arc_color = [255 255 191]/255;
line(real(arc1),imag(arc1),...
    'LineWidth',0.5,...
    'Color',arc_color);
line(real(arc2),imag(arc2),...
    'LineWidth',3,...
    'Color',arc_color);
 
axis equal

function [arc1,arc2] = triangleCurves(T)

arc1 = [];
arc2 = [];
for k = 1:height(T)
   t_k = T(k,:);
   switch t_k.Type
      case 'A'
         [arc1_k,arc2_k] = arcsA(t_k.Apex,t_k.Left,t_k.Right);
      case 'Ap'
         [arc1_k,arc2_k] = arcsAp(t_k.Apex,t_k.Left,t_k.Right);         
      case 'B'
         [arc1_k,arc2_k] = arcsB(t_k.Apex,t_k.Left,t_k.Right);
      case 'Bp'
         [arc1_k,arc2_k] = arcsBp(t_k.Apex,t_k.Left,t_k.Right);
   end
    arc1 = [arc1 arc1_k NaN];
    arc2 = [arc2 arc2_k NaN];
end

function [arc1,arc2] = arcsA(P1,P2,P3)
ray = P1 - P2;
theta1 = angle(ray);
theta2 = theta1 - deg2rad(72);
theta = linspace(theta1,theta2,10);
arc1 = P2 + 0.25 * abs(ray) * exp(1i * theta);

ray = P1 - P3;
theta1 = angle(ray);
theta2 = theta1 + deg2rad(72);
theta = linspace(theta1,theta2,10);
arc2 = P3 + 0.25 * abs(ray) * exp(1i * theta);

function [arc1,arc2] = arcsAp(P1,P2,P3)
ray = P1 - P2;
theta1 = angle(ray);
theta2 = theta1 - deg2rad(72);
theta = linspace(theta1,theta2,10);
arc2 = P2 + 0.25 * abs(ray) * exp(1i * theta);

ray = P1 - P3;
theta1 = angle(ray);
theta2 = theta1 + deg2rad(72);
theta = linspace(theta1,theta2,10);
arc1 = P3 + 0.25 * abs(ray) * exp(1i * theta);

function [arc1,arc2] = arcsB(P1,P2,P3)
ray = P1 - P2;
theta1 = angle(ray);
theta2 = theta1 - deg2rad(36);
theta = linspace(theta1,theta2,10);
arc2 = P2 + 0.75 * abs(ray) * exp(1i * theta);

ray = P1 - P3;
theta1 = angle(ray);
theta2 = theta1 + deg2rad(36);
theta = linspace(theta1,theta2,10);
arc1 = P3 + 0.25 * abs(ray) * exp(1i * theta);

function [arc1,arc2] = arcsBp(P1,P2,P3)
ray = P1 - P2;
theta1 = angle(ray);
theta2 = theta1 - deg2rad(36);
theta = linspace(theta1,theta2,10);
arc1 = P2 + 0.25 * abs(ray) * exp(1i * theta);

ray = P1 - P3;
theta1 = angle(ray);
theta2 = theta1 + deg2rad(36);
theta = linspace(theta1,theta2,10);
arc2 = P3 + 0.75 * abs(ray) * exp(1i * theta);
