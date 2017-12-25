%%
% Simple test for Frank-Wolfe algorithm


rep = '../results/frank-wolfe/';
[~,~] = mkdir(rep);

test = 2;

% square
K = [-1-1i; -1+1i; 1+1i; 1-1i];
% polygon
m = 3;
K = exp( 1i*pi/m + 2i*pi * (0:m-1)'/m );

% function 
switch test
    case 1
        p0 = 1+.3i; % target point
    case 2
        p0 = 0;
end

f = @(z)1/2*abs(z-p0).^2;
df = @(z)z-p0;

dotp = @(u,v)real(u*conj(v));

% draw
q = 256;
x = linspace(-1.05,1.05,q);
y = linspace(-.9,.9,q);
[Y,X] = meshgrid(y,x);
F = f(X+1i*Y);


niter = 18;
Z = [];

% initialize farthest possible
[~,i] = max(f(K));
z = K(i);
z = K(2);

clf; hold on;
imagesc(x,y,F');
plot(K([1:end 1]), 'k.-', 'MarkerSize', 25, 'LineWidth', 2);
for i=1:niter
    Z(i) = z;
    d = df(z);
    % find extremal point
    u = K(1); 
    for j=1:length(K)
        if dotp(d,K(j))<dotp(d,u)
            u = K(j);
        end
    end
    % update
    gamma = 2/(2+i);
    z = z + gamma*(u-z);
    % draw
    t = (i-1)/(niter-1);
    plot([Z(end) z], '.-', 'color', [t 0 1-t], 'LineWidth', 2, 'MarkerSize', 25);  
end

clf; hold on;
% background function
r = 20;
imagesc(x,y,F');
contour(x,y,F',r, 'k');
colormap(parula(r+1));
%
plot(K([1:end 1]), 'k.-', 'MarkerSize', 25, 'LineWidth', 2);
plot(real(p0), imag(p0), 'r.',  'MarkerSize', 25);
for i=1:niter-1
    t = (i-1)/(niter-2);
    plot(Z(i:i+1), '.-', 'color', [t 0 1-t], 'LineWidth', 2, 'MarkerSize', 25);    
end
axis off; axis equal
saveas(gcf, [rep 'fw-' num2str(test) '.png']);