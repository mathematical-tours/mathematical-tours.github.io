%%
% test for qr iterations

name = 'lapl';
name = 'filter';
name = 'randsym';
name = 'randn';

addpath('../toolbox/');
rep = MkResRep(name);

A = rand(4);

n = 200;
n = 30; 

t = 1:n; 
t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
A0 = X.^3-Y.^2;

% n = 4;

switch name
    case 'lapl'
        A0 = -2*diag(ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1); 
    case 'randn'
        A0 = randn(n);
    case 'randsym'
        A0 = randn(n);
        A0 = (A0+A0')/2;   
    case 'randeig'
        U = randn(n);
        A0 = U*diag(randn(n,1))*inv(U);
    case 'filter'
        A0 = 1./(.1+abs(X-Y));
end

% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));

niter = 300;
q = 60;
ndisp = unique(round(1+linspace(0,1,q).^3*(niter-1)));
k = 1;

A = A0; Z = eye(n);
for i=1:niter
    [Q,R] = qr(A);
    A = R*Q;
    [Z,Z1] = deal(Z*Q,Z);
    Z = Z * diag(sign(diag(Z*Z1'))) ;
    
    % Z(:,end) = Z(:,end) * sign(sum(Z(:,end)));
    % imagesc(Z(:,1:end-1)); axis image; axis off; drawnow;
    
    % apply colormap to the image
    if ndisp(k)==i
        U = log(1 + abs(A));
        % U = A;
        % imagesc(U); axis image; axis off; drawnow;
        imwrite(Upsc(rescale(U)*255,8), parula(256), [rep 'qr-' znum2str(k,2) '.png'] );
        imwrite(Upsc(rescale(Z(:,1:end-1))*255,8), parula(256), [rep 'eigv-' znum2str(k,2) '.png'] );
        k=k+1;
    end
    
    clf; imagesc(Q*R); drawnow;
    
end
% V = rescale(B);
% write_video(V, 'qr-movie', 'gif');