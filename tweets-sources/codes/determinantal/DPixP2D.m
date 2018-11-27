function X = DPixP2D(M,N,FC)
% [N,M] : image domain; FC : Fourier coefficients of the kernel C

%%%% Step 1 : Bernoulli sampling %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi = find(FC>rand(N,M));% Random sampling of the active frequencies
n = size(xi,1);  % Number of points of the DPP
xi2  = floor((xi-1)/N) + 1; % Actives frequencies column indices
xi1 = xi - (xi2 - 1) * N; % Actives frequencies row indices

str = sprintf(' Sample of size %d of a Determinantal Pixel Process', n); %titre de la figure

% Active eigenfunctions matrix
[X2, X1] = meshgrid(1:M,1:N);
X1 = reshape(X1,1,N*M);
X2 = reshape(X2,1,N*M);
V = exp(2*1i*pi*( (xi1*X1)/N + (xi2*X2)/M ) )/ sqrt(N*M);

%%%% Step 2 : Sampling the defined projection DPP %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(1,n);
e = zeros(n,n);

% Initialization
D = (1/n)*diag(V'*V);
X(1) = sum(cumsum(D) < rand(1))+1;

e(:,1) = V(:,X(1)) / sqrt(sum( abs(V(:,X(1))).^2 ));

% Sampling of the other points
for k = 2:n
    D = (diag(V'*V)' - sum(abs(e'*V).^2))/(n-k+1);
    X(k) = sum(cumsum(D) < rand(1))+1;
    w = V(:,X(k)) - sum(transpose((e'*V(:,X(k)))*ones(1,n)).*e,2);
    e(:,k) = w / sqrt(sum(abs(w).^2));
    
    
    % colormap
    DisplayDensity(reshape(D,[N M]), (k-2)/(n-2), X(1:k-1));
end

%%%% Vizualization of the sample %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Y = 255*ones(N,M); Y(X) = 0; 
% figure; imagesc(Y);colormap(gray);axis('square'); axis([0 N 0 M]); title(str);

end
