function [v,umax] = convert_distance_color(D,M, r,vmax)

% convert_distance_color - convert a distance function to a color image
%
%   A = convert_distance_color(D,M);
%
%   M is optional: background image.
%
%   Very useful to save a result of distance computation to an image file
%   with nice colors.
%
%   Copyright (c) 2007 Gabriel Peyre


n = size(D,1);
if nargin<2
    M = ones(n);
end

if nargin<3
r = 256;
r = 25;
end

c = parula(r);

U = D; 


umax= 1;
u = sort(U(U~=Inf));
if isempty(u)
    v = M; 
    return;
end
umax = u(round(end));
% umax = u(round(end));
if umax==0
    umax=1;
end

U(U==Inf) = min(U(U~=Inf));

% V = rescale(U); 
V = U/umax; V(V~=Inf) = min(V(V~=Inf),1);

% V = min(U/vmax,1);
I = floor((r-1)*V)+1;
v = c(I(:),:); 
v = reshape(v, [n n 3]);
A = repmat(rescale(M), [1 1 3]);
B = repmat(D, [1 1 3]);
v(B==Inf) = A(B==Inf);