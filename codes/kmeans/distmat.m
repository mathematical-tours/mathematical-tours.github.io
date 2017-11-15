function d = distmat(x,y, N0)

% distmat - compute pairwise L^2 distance matrix
%
%   d = distmat(x,y, N0);
%
%   d_ij = |x(:,i)-y(:,j)|
%
%   N0 is optional, and is used to sliced the computation by group of size
%       N0 to avoid memory overflow.
%
%   Copyright (c) 2015 Gabriel Peyre

if nargin<2
    y = x;
end


Nx = size(x,2);
Ny = size(y,2);

if nargin<3
    N0 = 5000;
end

if Nx>N0
    d = zeros(Nx,Ny);
    for k=1:ceil(Nx/N0)
        s = (k-1)*N0+1:min(Nx,k*N0);
        d(s,:) = distmat(x(:,s),y);
    end
    return;
end

d = sqrt( max(bsxfun(@plus,dot(x,x,1)',dot(y,y,1))-2*(x'*y), 0) );

%SLOWER
% a = sum(x.^2); b = sum(y.^2);
% d = repmat(b,Nx,1) + repmat(a',1,Ny) - 2*x'*y;
% d = sqrt(max(d,0));
        
end