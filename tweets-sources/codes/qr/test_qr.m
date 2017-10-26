%%
% test for qr iterations

A = rand(4);

N = 200;
U = randn(N);
A = U*diag(randn(N,1))*inv(U);

niter = 200;
B = [];
for i=1:niter
    [Q,R] = qr(A); 
    A = R*Q;
    U = log(1e-12 + abs(A));
    imagesc(U); axis image; axis off; drawnow;
   % if mod(i,2)==1
        B(:,:,end+1)=U;
   % end
end
V = rescale(B);
write_video(V, 'qr-movie', 'gif');