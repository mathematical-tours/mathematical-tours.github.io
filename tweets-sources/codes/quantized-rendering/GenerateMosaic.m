function f1 = GenerateMosaic(f,Z,kx)

n = size(Z,1);
m = size(Z,4);

% means
mu = mean(mean(Z,1),2);
mu = squeeze(mu)';

%
w = floor(size(f,1)/kx);
ky = floor(size(f,2)/w);
f = f(1:w*kx,1:w*ky,:,:);
k = kx*ky;
g = reshape(f, [w kx w ky 3]);
g = squeeze( mean(mean(g,1),3) );
G = reshape(g, [k 3]);

% randomization factor for NN search
q = 600;
% NN search
warning off;
I = knnsearch(G,mu,q);
warning on;

% randomize a little bit to avoid repetition
ind = floor(rand(k,1)*q)+1;
J = [];
for i=1:k
    J(i) = I(i,ind(i));
end
J = reshape(J,[kx ky]);

% re-create the final image
f1 = zeros(kx*n, ky*n,3);
for i=1:kx
    % progressbar(i,kx);
    for j=1:ky
        q = J(i,j); % selected patch
        x = (i-1)*n+1:i*n;
        y = (j-1)*n+1:j*n;
        a = double( Z(:,:,:,q) );
        a = permute(a, [2 1 3]);
        a = a + repmat( reshape(-mu(q,:),[1 1 3]) + g(i,j,:) , [n n] );
        a = clamp(a,0,255);
        f1(x,y,:) = a;
    end
end

end