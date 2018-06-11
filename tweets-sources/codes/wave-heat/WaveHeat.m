%%
% Wave vs. Heat equation.


addpath('../toolbox/');
rep = MkResRep();
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% step size
n = 256;


% Laplacian
a = [2:n n]; b = [1 1:n-1];
a = [2:n 1]; b = [n 1:n-1]; % per
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

% initialization
f0 = zeros(n);
f0(end/2,end/2) = 1;

t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

Z = [.7+.5i; .2+.8i; .4+.2i];
w = [1 -.6 .8];

s = .015;
f0 = zeros(n);
for k=1:length(Z)
    z = Z(k);
    f0 = f0 + w(k) * exp(-( (X-real(z)).^2+(Y-imag(z)).^2)/(2*s^2));
end

name = 'wave';
name = 'heat';
name = 'biheat';
% heat, df/dt = Delta(f)

switch name
    case 'heat'
        % heat, df/dt = Delta(f)
        niter = 7000;
        tau = .8;
    case 'biheat'
        % df/dt = Delta^2(f)
        niter = 7000;
        tau = .002;
    case 'wave'
        % heat, df^2/dt = Delta(f)
        niter = 7000;
        tau = .002;        
end



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
        saveas(gcf, [rep name '-' znum2str(k,2) '.png'], 'png');
    end
    switch name
        case 'heat'
           	f = f + tau * Delta(f);
        case 'biheat'
           	f = f + tau * Delta(Delta(f));
        case 'wave'
            % (2*f-f_new-f1) = tau*Delta(f)
            [f,f1] = deal( 2*f-f1+tau*Delta(f), f );
        case 'biwave'
            [f,f1] = deal( 2*f-f1+tau*Delta(Delta(f)), f );
    end
end
axis tight;


% AutoCrop(rep, [name '-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif