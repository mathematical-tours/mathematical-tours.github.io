% solve using Chambolle/Pock
% min_{p in S}?max_{q in S}?<Cp,q>
% min_{x in S}?max_{y in S}?<y,Kx> - f^*(y)+g(x)

name = 'smooth';
name = 'saddle';

addpath('../toolbox/');
rep = MkResRep(name);


% C-P
% x1 = x
% x <- Prox_{tau*g}(x - tau*K'*y)
% y <- Prox_{sigma*f^*}(y + sigma*K(2x-x1))

% gaussian blur
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

Proj = @(u)perform_simplex_projection(u',1)';

n = 50; p = 40;

n = 200; p = n; 
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

switch name
    case 'rand'        
        % random
        K = rand(p,n);
    case 'saddle'
        % unique saddle point
        K = -(X-.2).^2 + (Y-.6).^2;
    case 'smooth'
        % smooth field
        K = GFilt(randn(n),16/2);
        K = rescale(K,-1, 1);
end

K = K/max(abs(K(:)));


% display K
r = 15; % #levellines
clf; hold on;
imagesc(t,t,K');
contour(t,t,K',linspace(-1,1,r), 'k');
colormap(parula(r-1));
caxis([-1 1]);
axis image; axis off;
saveas(gcf, [rep 'K.png']);

niter = 5000;
q = 50; 
ndisp = round(1+(niter-1)*linspace(0,1,q).^3); 
ndisp = unique(ndisp);
k = 1;


tau = .05; sigma = .05; 
E = [];
x = ones(n,1)/n; % init
y = ones(p,1)/p;
for i=1:niter
    progressbar(i,niter);
    x1 = x;
    x = Proj(x - tau*K'*y);
    y = Proj(y + sigma*K*(2*x-x1));
    E(i) = norm(x-x1,1);
    if i==ndisp(k)
        s = (i-1)/(niter-1);
        clf;
        area(t,x, 'EdgeColor', [s 0 1-s], 'FaceColor', [s 0 1-s]); axis tight;
        set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/3 1]);
        saveas(gcf, [rep 'x-' znum2str(k,2) '.png']);
        %
        clf;
        area(t,y, 'EdgeColor', [s 0 1-s], 'FaceColor', [s 0 1-s]); axis tight;
        set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/3 1]);
        saveas(gcf, [rep 'y-' znum2str(k,2) '.png']);
        drawnow;
        k = k+1;
    end
end
clf;
plot(log10(E), '-');

% AutoCrop(rep, 'x-'); AutoCrop(rep, 'y-');

% via fictious play