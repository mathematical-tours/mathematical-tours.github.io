%%
% test for qr iterations

addpath('../toolbox/');
rep = MkResRep();

A = rand(4);

N = 200;
N = 30; 

U = randn(N);
A = U*diag(randn(N,1))*inv(U);

% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));


niter = 200;
B = [];
for i=1:niter
    [Q,R] = qr(A); 
    A = R*Q;
    U = log(1e-12 + abs(A));
    imagesc(U); axis image; axis off; drawnow;
    % apply colormap to the image
    if mod(i,2)==1
        imwrite(Upsc(rescale(U)*255,8), parula(256), [rep 'qr-' znum2str((i-1)/2+1,3) '.png'] );
    end
   % if mod(i,2)==1
	% B(:,:,end+1)=U;
   % end
end
% V = rescale(B);
% write_video(V, 'qr-movie', 'gif');