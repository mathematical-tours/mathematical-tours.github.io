%%
% Spacially varying wave vs. Heat equation.




name = 'iso';
name = 'anis';
name = 'ortho';
name = 'halfiso';

addpath('../toolbox/');
rep = MkResRep(name);


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% spatial step size
n = 256;

% Laplacian
a = [2:n n]; b = [1 1:n-1];
a = [2:n 1]; b = [n 1:n-1]; % per

%centered
Grad = @(f)cat(3, f(a,:)-f(b,:), f(:,a)-f(:,b))/2;
Div = @(v)( v(b,:,1)-v(a,:,1) + v(:,b,2)-v(:,a,2) )/2;


% grad
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f)/2;
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) )/2;



Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

f = randn(n);
mynorm = @(x)norm(x(:));
% should be 0
mynorm(-Div(Grad(f))-Delta(f))/mynorm(f);

% initialization
t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);

% orthogonal system
D = sqrt(X.^2+Y.^2)+1e-10;
V = cat(3,X./D,Y./D);
W = cat(3,-Y./D,X./D);

stensor = @(U)cat(3, U(:,:,1).^2, U(:,:,2).^2, U(:,:,1).*U(:,:,2));
tens_mult = @(T,v)cat(3, T(:,:,1).*v(:,:,1)+T(:,:,3).*v(:,:,2), T(:,:,3).*v(:,:,1)+T(:,:,2).*v(:,:,2) );


alpha = [.3 1]; % aniso

switch name
    case 'iso'
        alpha = [1 1];
    case 'anis'
        alpha = [.05 1];
    case 'halfiso'
        alpha = [.3 1];
    case 'ortho'
        alpha = [1 .01];
end

T = alpha(1)*stensor(V) + alpha(2)*stensor(W);

options.sub = round(n/20);
clf;
plot_tensor_field(T, [], options);
axis equal; axis tight; axis off;
saveas(gcf, [rep 'tensors-' name '.png'], 'png');
return;

Tmax = 150; 
tau = .1;
niter = round(Tmax/tau);

q = 50; 
disp_list = round( linspace(1,niter,q) );
kdisp = 1;

f  = randn(n);
for i=1:niter    
    if i==disp_list(kdisp)
        % display 
        s = std(f(:));
        r = 1;
        g = clamp( f/s, -r, r )/r;
        %
        r = 16; % #levellines
        clf; hold on;
        u = linspace(0,1,n);
        imagesc(u,u,g);
        % contour(u,u,g,linspace(-1,1,r), 'k');
        colormap(parula(r-1));
        caxis([-1 1]);
%        contour(u,u,Rho,[1 1]*mean(Rho(:)), 'b--', 'LineWidth', 2);        
        axis image; axis off; axis ij;
        drawnow;
        %
        kdisp = kdisp+1;
        % saveas(gcf, [rep name '-' znum2str(kdisp,2) '.png'], 'png');
        imwrite((g+1)/2, [rep name '-' znum2str(kdisp,2) '.png'], 'png');
    end
    f = f - tau * Div( tens_mult(T,Grad(f)) );
end
axis tight;


% AutoCrop(rep, [name '-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif