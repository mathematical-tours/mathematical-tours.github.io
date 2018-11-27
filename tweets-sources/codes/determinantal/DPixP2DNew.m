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
V = exp(2i*pi*( (xi1*X1)/N + (xi2*X2)/M ) )/ sqrt(N*M);

%%%% Step 2 : Sampling the defined projection DPP %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(1,n);
e = zeros(n,n);

% Initialization
if 0
    D = (1/n)*diag(V'*V);
else
    D = zeros(N*M,1);
    for k=1:n
        D = D +  abs(V(k,:)').^2;
    end
    D = D/n;
end

X(1) = sum(cumsum(D) < rand(1))+1;

e(:,1) = V(:,X(1)) / sqrt(sum( abs(V(:,X(1))).^2 ));

% Sampling of the other points
for k = 2:n
    if 0
        D = (diag(V'*V)' - sum(abs(e'*V).^2))/(n-k+1);
    else
%         D = zeros(N*M,1);
%         for r=1:n
%             D = D +  abs(V(r,:)').^2;
%         end
        D = n/(N*M);
        D = ( D - sum(abs(e'*V).^2)' ) /(n-k+1);
    end
    %% Display %%    
    % colormap
    DisplayDensity(reshape(D,[N M]), (k-2)/(n-2), X(1:k-1));
    % D = [ones(N,N/2), zeros(N,N/2)]; D = D(:)/sum(D(:));

    %
    X(k) = sum(cumsum(D) < rand(1))+1;
    w = V(:,X(k)) - sum(transpose((e'*V(:,X(k)))*ones(1,n)).*e,2);
    e(:,k) = w / sqrt(sum(abs(w).^2));
end

end
