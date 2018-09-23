function [g,nj] = PerformHaarSphere(f,direction,J,nj)



% Precompute the local wavelet matrix, which contains the local vector and
% three orthognal detail directions.
randn('seed', 1);
U = randn(4);
U(:,1) = 1;
[U,R] = qr(U);

if direction==+1
    % FWD
    g = f;
    nj = length(f);
    for j=1:J-1
        fj = g(1:nj);
        fj = reshape(fj, [nj/4 4]);
        fj = fj*U;
        g(1:nj) = fj(:);
        nj = nj/4;
    end
else
    % BWD
    g = f;
    for j=1:J-1
        fj = g(1:4*nj);
        fj = reshape(fj, [nj 4]);
        fj = fj*U';
        g(1:4*nj) = fj(:);
        nj = nj*4;
    end
end