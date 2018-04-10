%%
% Wave vs. Heat equation.


addpath('../toolbox/');
rep = MkResRep();
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

s = .01;
f0 = zeros(n);
for k=1:length(Z)
    z = Z(k);
    f0 = f0 + w(k) * exp(-( (X-real(z)).^2+(Y-imag(z)).^2)/(2*s^2));
end


niter = 1000;
tau = .1;

H = 4;


% svg run
ndisp = max(1,ceil(niter/50)); k = 0; 
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
        contour(u,u,g,linspace(-1,1,r), 'k');
        colormap(jet(r-1));
        caxis([-1 1]);
        axis image; axis off; axis ij;
        %
        saveas(gcf, [rep  'porous-' num2str(H) '-' znum2str(k,2) '.png'], 'png');
    end
    f = f + tau * Delta(f).^H;
end
axis tight;


% AutoCrop(rep, [name '-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif