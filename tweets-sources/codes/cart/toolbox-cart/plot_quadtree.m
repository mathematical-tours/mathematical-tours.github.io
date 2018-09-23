function plot_quadtree(W, f, options)

% plot_quadtree - plot an image quadtree
%
%   plot_quadtree(T, f, options);
%
%   f is a background image.
%
%   Copyright (c) 2010 Gabriel Peyre


options.null = 0;
if nargin<2
    f = [];
end

n = size(f,1);
J = length(W);
str = 'r';
str_geom = 'b';


hold on;

% display image
if ~isempty(f)
    %f = f';
    %f = f(end:-1:1,:);
    %f = fliplr(f);
    imagesc([0 1],[0 1],f);
    colormap gray(256);
end
plot_square([0,0], 1, str);
axis square;
axis off;
axis equal;

cx = [.5];
cy = [.5];
v = 1;
for j=1:J
    z = v(:)*0; z(v(:)) = 1:length(v(:));
    w = 1/2^j; % width of a square
    for k=1:length(W{j})
        if W{j}(k)==0
           plot_cross([cy(z(k)) cx(z(k))], w*2, str);
        end
    end
    % update for the next scale
    cx1 = zeros(2^j,2^j);
    cx1(1:2:end,1:2:end) = cx-w/2;
    cx1(2:2:end,1:2:end) = cx+w/2;
    cx1(1:2:end,2:2:end) = cx-w/2;
    cx1(2:2:end,2:2:end) = cx+w/2; 
    cy1 = zeros(2^j,2^j); 
    cy1(1:2:end,1:2:end) = cy-w/2;
    cy1(2:2:end,1:2:end) = cy-w/2;
    cy1(1:2:end,2:2:end) = cy+w/2;
    cy1(2:2:end,2:2:end) = cy+w/2;   
    cx = cx1;   
    cy = cy1;   
    % update for the next scale
    v1 = zeros(2^j,2^j);
    v1(1:2:end,1:2:end) = v*4-3;
    v1(2:2:end,1:2:end) = v*4-2;
    v1(1:2:end,2:2:end) = v*4-1;
    v1(2:2:end,2:2:end) = v*4;    
    v = v1;
end

axis ij

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_cross(pos, w, str)

pos = swap_pos(pos);

if nargin<3
    str = 'r';
end

x = [pos(1)-w/2, pos(1)+w/2];
y = [pos(2), pos(2)];
plot(x,y, str);
x = [pos(1), pos(1)];
y = [pos(2)-w/2, pos(2)+w/2];
plot(x,y, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_square(pos, w, str)

% pos = swap_pos(pos);

if nargin<3
    str = 'r';
end

x = [pos(1), pos(1)+w, pos(1)+w, pos(1), pos(1)];
y = [pos(2), pos(2), pos(2)+w, pos(2)+w, pos(2)];
plot(x,y, str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_square_geometry(theta,pos,w, str)

if nargin<4
    str = 'b';
end

% pos = pos(2:-1:1);
x = pos(1)+w/2 + w/2*[cos(theta), -cos(theta)];
y = pos(2)+w/2 + w/2*[sin(theta), -sin(theta)];
plot(x,y, str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos1 = swap_pos(pos)

pos1 = pos;
% pos1 = pos(2:-1:1);
%% pos1(1) = 1-pos1(1);
