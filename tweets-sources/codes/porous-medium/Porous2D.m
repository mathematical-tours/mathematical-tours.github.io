%%
% Non linear heat equation.


addpath('../toolbox/');
rep = MkResRep('2d');
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% step size
n = 256/2;

% Laplacian
a = [2:n n]; b = [1 1:n-1];
a = [2:n 1]; b = [n 1:n-1]; % per
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

% initialization
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
%
Z = [.36+.47i; .15+.3i];
w = [1;.7];

k = 10;
rand('state', 123);
Z = rand(k,1) + 1i*rand(k,1);
w = 1 + .2*rand(k,1);

s = .04;
f0 = zeros(n);
for k=1:length(Z)
    z = Z(k);
    f0 = f0 + w(k) * exp(-( (X-real(z)).^2+(Y-imag(z)).^2)/(2*s^2));
end



H = 3;
H = 1;
switch H
    case 4
        Tmax = 200000;
        tau = .9;
    case 3
        Tmax = 5000;
        tau = .2;
    case 2
        Tmax = 10000;
        tau = .9;
    case 1
        Tmax = 500;
        tau = .6;
end
niter = round(Tmax/tau);



% svg run
ndisp = max(1,ceil(niter/50)); 
k = 0; 
%
f  = f0;
f1 = f0;
for i=1:niter    
    if mod(i,ndisp)==1
        k = k+1;
        % display 
        g = f/max(abs(f(:)));
        r = 16; % #levellines
        clf; hold on;
        u = linspace(0,1,n);
        imagesc(u,u,g);
        contour(u,u,g,linspace(0,1,r), 'k');
        colormap(parula(r-1));
        caxis([0 1]);
        axis image; axis off; axis ij;
        drawnow;
        %
        saveas(gcf, [rep  'porous-' num2str(H) '-' znum2str(k,2) '.png'], 'png');
    end
    f = f + tau * Delta(f.^H);
end
axis tight;


% AutoCrop(rep, ['porous-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif