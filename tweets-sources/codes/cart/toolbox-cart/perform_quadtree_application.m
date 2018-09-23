function y = perform_quadtree_application(y,T,f)

% perform_quadtree_application - apply function to quadtree
%
%   y = perform_quadtree_application(y,T,f);
%
%   The function b=f(a) is applied to a matrix a of size (w,w,p) of p
%   patches of size (k,k), and should return a matrix of size (w,w,p).
%
%   Copyright (c) 2010 Gabriel Peyre


n = size(y,1);
J = length(T);

v = 1;
for j=1:J
    
    % cut squares of size n/2^(j-1) x n/2^(j-1)
    k = 2^(j-1);
    w = n/k;
    [dY,dX] = meshgrid(0:w-1,0:w-1);
    dX = repmat( reshape(dX,[1 1 w w]), [k k 1] );
    dY = repmat( reshape(dY,[1 1 w w]), [k k 1] );
    x = 1:w:n;
    [Y,X] = meshgrid(x,x);
    X = repmat(X, [1 1 w w]);
    Y = repmat(Y, [1 1 w w]);
    I = X+dX + (Y+dY-1)*n;
    % patches of size (k,k,w,w);
    P = y(I);
    P = permute(P, [3 4 1 2]);
    P = reshape(P, [w w k*k]);     
    % leaf patches    
    J = find( T{j}(v(:))==1 );    
    P(:,:,J) = f(P(:,:,J));
    P = reshape(P, [w w k k]);
    P = permute(P, [3 4 1 2]);
    % 
    y(I) = P;
    % update for the next scale
    v1 = zeros(2^j,2^j);
    v1(1:2:end,1:2:end) = v*4-3;
    v1(2:2:end,1:2:end) = v*4-2;
    v1(1:2:end,2:2:end) = v*4-1;
    v1(2:2:end,2:2:end) = v*4;    
    v = v1;
end

return;


% adaptive quadtree thresholding
J = length(T);
cx = [1];
cy = [1];
v = 1;
for j=1:J
    z = v(:)*0; z(v(:)) = 1:length(v(:));
    % w = 1/2^j; % width of a square
    w = n/2^(j-1);
    for k=1:length(T{j})
        if T{j}(k)==1
            % position of the center
            selx = cx(z(k)):cx(z(k))+w-1;
            sely = cy(z(k)):cy(z(k))+w-1;
            y(selx,sely) = f(y(selx,sely));
        end
    end
    % update for the next scale
    cx1 = zeros(2^j,2^j);
    cx1(1:2:end,1:2:end) = cx;
    cx1(2:2:end,1:2:end) = cx+w/2;
    cx1(1:2:end,2:2:end) = cx;
    cx1(2:2:end,2:2:end) = cx+w/2;
    cy1 = zeros(2^j,2^j);
    cy1(1:2:end,1:2:end) = cy+w/2;
    cy1(2:2:end,1:2:end) = cy+w/2;
    cy1(1:2:end,2:2:end) = cy;
    cy1(2:2:end,2:2:end) = cy;
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