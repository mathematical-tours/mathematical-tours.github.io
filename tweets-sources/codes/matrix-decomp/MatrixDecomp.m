
addpath('../toolbox/');
rep = MkResRep();

% example of matrix decomposition.
n = 20;


% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));
% turn gray scale image into color image
CM = parula(256);
TurnColor = @(C)reshape(CM(1+floor(rescale(Upsc(C,8))*255),:), [8*n 8*n 3]);
%
mysave = @(A,str,it)imwrite(TurnColor(A), [rep str '-' znum2str(it,2) '.png']);

q = 50;

vlist = .05 + 100*linspace(0,1).^3; 

for it=1:q

[Y,X] = meshgrid(1:n,1:n);
A = 1./(abs(X-Y)+vlist(it));

% A = P'*L*U;
[L,U,P] = lu(A);

% A=Q*R
[Q,R] = qr(A);

% A=C'*C
C = chol(A);

clf; 
subplot(2,3,1); imagesc(A); title('A');
subplot(2,3,2); imagesc(L); title('L');
subplot(2,3,3); imagesc(U); title('U');
subplot(2,3,4); imagesc(Q); title('Q');
subplot(2,3,5); imagesc(R); title('R');
subplot(2,3,6); imagesc(C); title('C');
drawnow;

mysave(A,'A',it);
mysave(L,'L',it);
mysave(U,'U',it);
mysave(Q,'Q',it);
mysave(R,'R',it);
mysave(C,'C',it);


end