%%
% Helmholtz decomposition 

addpath('../toolbox/');
rep = MkResRep();

n = 256;

% gaussian blur
t = [0:n/2,-n/2+1:-1]';
G = @(s)exp(-t.^2/(2*s^2)); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));



% grad and div 
a = [2:n n]; b = [1 1:n-1]; % symmetric bc
a = [2:n 1]; b = [n 1:n-1]; % periodic bc
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f);
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) );
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) );

% should be 0
f = randn(n);
mynorm = @(x)norm(x(:));
mynorm(-Div(Grad(f))-Delta(f))/mynorm(f)

% Fourier transform of Delta
[Y,X] = meshgrid(0:n-1,0:n-1);
mu = sin(X*pi()/n).^2; mu = -4*( mu+mu' );
% check correctness
f = randn(n);
mynorm( fft2(f).*mu - fft2(Delta(f))) / mynorm(f)

% pseudo-inverse of Laplacian
mu(1) = 1; % avoid 0 division
DeltaInv = @(u)real( ifft2( fft2(u) ./ mu ) );
% Projection on incompressible flows.
ProjI = @(u)u + Grad(DeltaInv(Div(u)));


% generate random gaussian fields

s = 10;
u = cat(3, GFilt(randn(n),s), GFilt(randn(n),s) );
v = ProjI(u);
w = u-v;

mynorm(Div(v))/mynorm(v)
mynorm(Curl(w))/mynorm(w)


sf = 1.2; % stretch factor
lw = 2;
q = 20;
k = round(n/q); % subsampling
plotvf = @(v,col)quiver(v(1:k:end,1:k:end,2),v(1:k:end,1:k:end,1), sf, col, 'LineWidth', lw);
%
V = {u v w};
Col = {'k' 'b' 'r'};
names = {'original' 'divfree' 'curlfree'};
for i=1:length(names)
    clf;
    plotvf(V{i}, Col{i});
    axis equal; axis off;
    saveas(gcf, [rep names{i} '-' num2str(s) '.eps'], 'epsc');
end

%% 
% Generate videos

repv = [rep 'video/'];
[~,~] = mkdir(repv);

q = 50; % #frames
slist = linspace(2,150,q);
u0 = randn(n,n,2);
for is=1:q
    s = slist(is);
    u = cat(3, GFilt(u0(:,:,1),s), GFilt(u0(:,:,2),s) );
    v = ProjI(u); w = u-v;
    V = {u v w};
    for i=1:length(names)
        clf;
        plotvf(V{i}, Col{i});
        axis equal; axis off;
        saveas(gcf, [repv names{i} '-' znum2str(is,2) '.png'], 'png');
    end
end


% create gifs
% AutoCrop(repv, ['original-1-'])
% Use % convert interp-*.png interp.gif to generate the gif using imagemagik


