%%
% Test for projected gradient descent.

name = 'square';
name = 'polygon';
name = 'disk';

addpath('../toolbox/');
rep = MkResRep(name);

n = 256; % for display

t = linspace(-2,2,n);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;

% constraint
switch name
    case 'polygon'
        p = 5;
        E = exp(2i*pi*(0:p-1)/p);
        s = linspace(0,1,100);
        Cc = [];
        for i=1:length(E)
            j = mod(i,length(E))+1;
            Cc = [Cc, (1-s)*E(i)+s*E(j)];
        end
        C = reshape( inpolygon(Y(:),X(:),real(E),imag(E)), [n n]);         
        Proj = @(z)ProjPoly(z,Cc,E);
    case 'disk'
        C = abs(Z)<=1;
        Proj = @(z)(abs(z)<1).*z + (abs(z)>=1) .* (z./abs(z));
        Cc = exp(2i*pi*linspace(0,1,200));
        % Extremal points
        E = Cc;
    case 'square'
        rho = .8;
        linf = @(z)max(abs(real(z)),abs(imag(z)));
        C = linf(Z)<=rho;
        Proj = @(z)(linf(z)<rho).*z + rho*(linf(z)>=rho) .* (z./linf(z));
        Cc = exp(2i*pi*linspace(0,1,200));
        Cc = rho*Cc./linf(Cc);
        E = rho*[1+1i,1-1i,-1+1i,-1-1i];
    case 'diamond'
        rho = 1.5;
        lun = @(z)abs(real(z)) + abs(imag(z));
        C = lun(Z)<=rho;
        Proj = @(z)(lun(z)<rho).*z + rho*(lun(z)>=rho) .* (z./lun(z));
        Cc = exp(2i*pi*linspace(0,1,200));
        Cc = rho*Cc./lun(Cc);
        E = rho*[1,-1,1i,-1i];
end

% function
z0 = 0+1.4i;
z0 = 1.3+1.4i;
%
F = 1/2*abs(Z-z0).^2;
F = sqrt(F/max(F(:)));
Grad = @(z)z-z0;
%
m = [4 1];
F = 1/2*( m(1)*real(Z-z0).^2 + m(2)*imag(Z-z0).^2 );
Grad = @(z)m(1)*real(z-z0) + m(2)*1i*imag(z-z0);
Grad = @(z)z-z0;
%
tau = .4;

% render the constraint;
r = 12;
F = sqrt(F/max(F(:)));
CM = parula(r-1);
I = min( 1+floor(rescale(F')*(r-1)), r-1);
U = reshape(CM(I,:), [n n 3]);
J = find(C==1);
U([J,J+n*n,J+2*n*n]) = U([J,J+n*n,J+2*n*n])/2;

q = 50;
for it=1:q
    
    s = (it-1)/(q-1);
    z0 = 1.7*exp(2i*pi*s);
    % run projected descent
    niter = 30; z = z0;
    for i=1:niter
        z(end+1) = z(end)-tau*Grad(z(end));
        z(end+1) = Proj(z(end));
    end
    % run frankwolf
    niter = 60; w = z0;
    for i=1:niter
        [~,j] = min(real(conj(Grad(w(end)))*E));
        gamma = 2/(2+i-1);
        w(end+1) = w(end) + gamma*(E(j)-w(end));
    end
    
    
    % display
    clf; hold on;
    imagesc(t,t,U);
    contour(t,t,F',linspace(0,1,r), 'k');
    colormap(parula(r-1));
    caxis([0 1]);
    plot(Cc,'k', 'LineWidth', 2);
    plot(z, 'r.-', 'LineWidth', 2, 'MarkerSize', 25);
    plot(w, 'g.-', 'LineWidth', 2, 'MarkerSize', 25);
    axis image; axis off;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end

% AutoCrop(rep, 'anim');


