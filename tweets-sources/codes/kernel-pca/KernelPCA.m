%%
% Schölkopf, Bernhard (1998). "Nonlinear Component Analysis as a Kernel Eigenvalue Problem"

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

% generate point in anulus
ann = @(r,n) ( r + .05*randn(n,1) ) .* exp(2i*pi*(rand(n,1)));

n = 1400;
n = 700;
x = ( .6 + .1*randn(n,1) ) .* exp(.7*2i*pi*(+5 + rand(n,1)));





% two anulus
s = .6;
n1 = round(n*s); n2 = n-n1;
x = [ann(.8,n1); ann(.0,n2)];
I = 1:n1; J = n1+1:n;

% three anulus
s = .33;
n1 = round(n*s); n2 = n1; n3 = n-(n1+n2);
x = [ann(.8,n1); ann(.4,n2); ann(.0,n3)];
I = 1:n1; J = n1+1:n;





% two moons
n1 = n/2;
x0 = ann(1,n);
x = x0;
m = median(imag(x));
I = find(imag(x)<m); J = find(imag(x)>=m);
x(I) = x(I)+1 +.35i;
x = rescale(real(x),-1,1) + rescale(imag(x),-.7,.7)*1i;
x = x*.8;


clf; hold on;
plot(x(I), 'b.'); plot(x(J), 'r.');

% centering matrix
Jc = eye(n) - ones(n)/n;
gaussian = @(x,y,s)exp(-abs( x-transpose(y) ).^2/(2*s^2));

% kernel
s = .2;


% grid of point
p = 150;
t = linspace(-1,1,p);
[Y,X] = meshgrid(t,t);
z = X(:)+1i*Y(:);

q = 50;
slist = linspace(.05,.5,q);
Uold = randn(n);
for it=1:q
    s = slist(it);
    K = gaussian(x,x,s);
    [U,D] = eig(Jc*K*Jc); D = diag(real(D));
    U = U * diag(sign(diag(U'*Uold))); % realign
    
    % project point
    x1 = K*( U(:,1) +  1i*U(:,2) );
    x1 = x1/std(abs(x1));
    clf; hold on;
    plot(x1(I), '.', 'color', [0 .8 0]); plot(x1(J), 'r.');
    axis([-1 1 -1 1]*5); axis off;
    axis tight; 
    drawnow;
    Uold = U;   
    mysaveas('proj',it);
end

save_mode = 1;

Uold = randn(n);
r = 15; % #levellines
for it=1:q
    s = slist(it);
    
    K = gaussian(x,x,s);
    [U,D] = eig(Jc*K*Jc); D = diag(real(D));
    U = U * diag(sign(diag(U'*Uold))); % realign
    
    z1 = gaussian(z,x,s)* U(:,1:10);
    z1 = reshape(z1,p,p, 10);
    
    if not(save_mode)
        clf;
    end
    for k=1:4
        if not(save_mode)
            subplot(2,2,k);
        else
            clf;
        end
        hold on;
        R = z1(:,:,k);
        R = R/max(abs(R(:))); R = (R+1)/2;
        imagesc(t,t,R');
        contour(t,t,R',linspace(0,1,r), 'k');
        colormap(parula(r-1));
        caxis([0 1]);
        plot(x(I), '.', 'color', [0 .8 0]); plot(x(J), 'r.');
        axis image; axis off;
        mysaveas(['eig' num2str(k)],it);
    end
    drawnow;
    
end