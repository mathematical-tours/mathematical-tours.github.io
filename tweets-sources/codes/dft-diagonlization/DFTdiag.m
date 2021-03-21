%%
% Diagonalize the FFT
% TODO: compute the multiplicies
% https://arxiv.org/pdf/1501.07646.pdf

% Hermite functions ? the eigenfunctions of the continuous Fourier transform. 

%  B. W. Dickinson and K. Steiglitz.  Eigenvectors and functions of the discrete Fourier transform.
% F. A. Grunbaum. The eigenvectors of the discrete Fourier transform

% F*F' = F'*F = Id

% F'=F^{-1}
n = 32;
[Y,X] = meshgrid(0:n-1,0:n-1);
F = 1/sqrt(n) * exp(2i*pi/n * X.*Y);

% problem : here U' is not equal to inv(U) -- but it could be ... 
[U,D] = eig(F); D = diag(D);
norm(F - U*diag(D)*inv(U), 'fro') / norm(F,'fro')
[D,I] = sort(angle(D)); U = U(:,I);

% F = U*D*U' U*U'--> F*F = 

A = logm(F)./1i;