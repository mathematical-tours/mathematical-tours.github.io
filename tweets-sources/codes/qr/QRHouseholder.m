%%
% Display QR factorization progress


addpath('../toolbox/');
rep = MkResRep();

% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));

n = 12;
A0 = randn(n);

        
        
imwrite(rescale(Upsc(A0,4)), [rep 'A.png']);


% Householder reflections for QR decomposition.
% [R,U] = house_qr(A) returns

A = A0;

H = @(u,x) x - u*(u'*x);
[m,n] = size(A);
U = zeros(m,n);
R = A;
for j = 1:min(m,n)
    u = house_gen(R(j:m,j));
    U(j:m,j) = u;
    R(j:m,j:n) = H(u,R(j:m,j:n));
    R(j+1:m,j) = 0;
    imagesc(R);
    % A=QR
    Q = A*inv(R);
    % imagesc(Q);
    drawnow;
    imwrite(rescale(Upsc(Q,4)), [rep 'Q-' znum2str(j,2) '.png']);
    imwrite(rescale(Upsc(R,4)), [rep 'R-' znum2str(j,2) '.png']);
end

%% Same but using Givens rotations

A = A0;
[m,n] = size(A);
Q = eye(m);
R = A;

imwrite(rescale(Upsc(A0,8)), [rep 'A.png']);
it = 1;
for j = 1:n
    for i = m:-1:(j+1)
        G = eye(m);
        [c,s] = givensrotation( R(i-1,j),R(i,j) );
        G([i-1, i],[i-1, i]) = [c -s; s c];
        R = G'*R;
        Q = Q*G;
        imagesc(Q); drawnow;
        % imagesc(Q); drawnow;
        imwrite(rescale(Upsc(Q,8)), [rep 'Q-' znum2str(it,3) '.png']);
        imwrite(rescale(Upsc(R,8)), [rep 'R-' znum2str(it,3) '.png']);
        it = it+1;
    end
end

