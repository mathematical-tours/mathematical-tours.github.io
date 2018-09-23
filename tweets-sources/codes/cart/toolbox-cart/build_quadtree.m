function T = build_quadtree(f, J, err)

% build_quadtree - build a quadtree penalty for an image
%
%   T = build_quadtree(f, J, err);
%
%   f is an (n,n) pixel image
%   J is the number of scales of the tree
%   err=@(P) is a callback that takes as input a 
%       set of patches P of size (k,k,w^2) and
%       output a set of penalization of size (k,k).
%
%   The indexing of the nodes of the tree reads:
%
% 1
%
% 1|2
% -+-
% 3|4
%
%  1 2 | 5 6
%  3 4 | 7 8
% -----+-----
% 9  10|13 14 
% 11 12|15 16
%
%   Copyright (c) 2010 	Gabriel Peyre
%			Jalal Fadili



n = size(f,1);
q = size(f,3);
T = {};
v = 1;
for j=1:J
    % cut squares of size n/2^(j-1) x n/2^(j-1)
    k = 2^(j-1);
    w = n/k;
    % [dY,dX] = meshgrid(0:w-1,0:w-1);
    [dX,dY,dZ] = ndgrid(0:w-1,0:w-1,1:q);
    dX = repmat( reshape(dX,[1 1 w^2*q]), [k k 1] );
    dY = repmat( reshape(dY,[1 1 w^2*q]), [k k 1] );
    dZ = repmat( reshape(dZ,[1 1 w^2*q]), [k k 1] );
    x = 1:w:n;
    [Y,X] = meshgrid(x,x);
    X = repmat(X, [1 1 w^2*q]);
    Y = repmat(Y, [1 1 w^2*q]);
    I = X+dX + (Y+dY-1)*n + (dZ-1)*n*n;
    % patches
    P = f(I);
    % compute error of approximation
    u = zeros(k,k,q);
    for iq=1:q
      u(:,:,iq) = err(P(:,:,(iq-1)*w*w+1:iq*w*w),w*w);
    end
    mu = mean(u,3);
    T{j}(v(:)) = mu;
    T{j} = T{j}(:);
    % update for the next scale
    v1 = zeros(2^j,2^j);
    v1(1:2:end,1:2:end) = v*4-3;
    v1(2:2:end,1:2:end) = v*4-2;
    v1(1:2:end,2:2:end) = v*4-1;
    v1(2:2:end,2:2:end) = v*4;    
    v = v1;
end
