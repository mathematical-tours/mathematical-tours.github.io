% NMF:
% V ~ WH
% min_{W>=0,H>=0, W(:,i) in simplex} KL(V|WH)

% PCA constrains the columns of W to be orthonormal and the rows of H to be orthogonal to each other. 

% W <- W . [ ( V/(W*H) ) * H^T ]
% W <- ProjRows(W)
% H <- H . [ W^T * ( V/(W*H) )  ]

suma = @(x)sum(x(:));
KL = @(V,Z)suma( V .* log(V./Z)-V+Z );

n = 20; % dimension of data
p = 100; % #data
s = 5; % #hidden factors
s = 15*3;

V = rand(n,p);
W = rand(n,s); 
H = rand(s,p);

niter=50; E = [];
for i=1:50
    W = W .* ( ( V./(W*H) ) * H' );
    W = W ./ sum(W,1); % simplex constraint
    H = H .* ( W' * ( V./(W*H) ) );
    E(i) = KL(V,W*H);
end
V1 = W*H;

clf;
plot(log10(E/E(1)));