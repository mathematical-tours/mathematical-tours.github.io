%%
% Spacially varying wave vs. Heat equation.


addpath('../toolbox/');
rep0 = MkResRep();
rep = [rep0 'aniso/'];
[~,~] = mkdir(rep);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% step size
n = 256;


% Laplacian
a = [2:n n]; b = [1 1:n-1];
a = [2:n 1]; b = [n 1:n-1]; % per

% grad
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f)/2;
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) )/2;
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

f = randn(n);
mynorm = @(x)norm(x(:));
% should be 0
mynorm(-Div(Grad(f))-Delta(f))/mynorm(f)

% initialization
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

Z = [.7+.5i; .2+.8i; .4+.2i];
w = [1 -.6 .8];

Z = [.5+.5i];
Z = [.24+.47i; .8+.73i];
%
Z = [.36+.47i; .41+.9i];
w = [1 -.4];
%
Z = [.36+.47i];
w = [1];
%
Z = [.36+.47i; .15+.3i];
w = [1;-.7];

s = .01;
f0 = zeros(n);
for k=1:length(Z)
    z = Z(k);
    f0 = f0 + w(k) * exp(-( (X-real(z)).^2+(Y-imag(z)).^2)/(2*s^2));
end

% weight 
d = .2;
Rho = .2 + ( (abs(X-.5)<d) & (abs(Y-.5)<d) ) ;

rho_name = 'room';
Rho = load_image(rho_name, n);
Rho = rescale(sum(Rho, 3))>.5;
if Rho(1)==1
    Rho = 1-Rho;
end
Rho = .0001+Rho;

DeltaRho = @(f)-Div(Grad(f) .* repmat(Rho, [1 1 2]) );

name = 'heat';
name = 'wave';

switch name
    case 'heat'
        % heat, df/dt = Delta(f)
        niter = 4*7000;
        tau = .95;
        fend = .2; 
    case 'wave'
        % heat, df^2/dt = Delta(f)
        niter = 7000;
        tau = .002; 
        fend = .1; finit = .3;       
end



% svg run
ndisp = max(1,ceil(niter/60)); k = 0; 
%
f  = f0;
f1 = f0;
for i=1:niter    
    if mod(i,ndisp)==1
        k = k+1;
        % display 
        r = (i-1)/(niter-1);
        g = f / ( (1-r)*finit + r*fend );
        g = f/max(abs(f(:)));
        %
        r = 16; % #levellines
        clf; hold on;
        u = linspace(0,1,n);
        imagesc(u,u,g);
        contour(u,u,g,linspace(-1,1,r), 'k');
        colormap(jet(r-1));
        caxis([-1 1]);
        contour(u,u,Rho,[1 1]*mean(Rho(:)), 'b--', 'LineWidth', 2);        
        axis image; axis off; axis ij;
        %
        saveas(gcf, [rep name '-' znum2str(k,2) '.png'], 'png');
    end
    switch name
        case 'heat'
           	f = f + tau * DeltaRho(f);
        case 'wave'
            % (2*f-f_new-f1) = tau*Delta(f)
            [f,f1] = deal( 2*f-f1+tau*DeltaRho(f), f );
    end
end
axis tight;


% AutoCrop(rep, [name '-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif