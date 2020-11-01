%%
% 2D networts


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 100;
s = [2:n,1]; t = [n,1:n-1];

% weights
T = linspace(-1,1,n);
[Y,X] = meshgrid(T,T);
Z = double( sqrt(X.^2+Y.^2)<=.7 );
Z = double( max(abs(X),abs(Y))<=.7 );

f = @(s)rescale(-abs(s), 0, 1);
W1 = f(Z-Z(s,:));
W2 = f(Z-Z(:,s));
W3 = f(Z-Z(t,:));
W4 = f(Z-Z(:,t));

W1 = W1*0+1;
W2 = W1*0+1;
W3 = W1*0+1;
W4 = W1*0+1;

% initial configuration 
Z = sign(randn(n));

niter = 200;
p = 3000; % #neurons to flip each time
E = []; 
j = 0;
for it=1:niter
    if 1 %mod(it,1)==1
    clf;
    imagesc(Z); 
    colormap gray;
    axis image; axis off;
    drawnow;
    j = j+1;
    mysaveas(j);
    end    
    %
    A = W1 .* Z(s,:) + W2 .* Z(:,s) + W3 .* Z(t,:) + W4 .* Z(:,t) + ...
         + 1e-10*randn(n); % + ...
       % X(s,s) + X(t,t) + X(s,t) + X(t,s);
    U= W1.*(Z-Z(s,:)).^2 + W2.*(Z-Z(:,s)).^2 + W3.*(Z-Z(t,:)).^2 + W4.*(Z-Z(:,t)).^2;
    E(it) = mean(U(:));
    I = randperm(n*n); I = I(1:p);
    Z(I) = sign(A(I)); 
end


