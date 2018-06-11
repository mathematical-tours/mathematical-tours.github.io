% compare heat to tv flow

eqname = 'heat';
eqname = 'tv';

addpath('../toolbox/');
rep = MkResRep(eqname);

n = 256;
name = 'hibiscus';
f0 = rescale(sum(load_image(name, n),3));
[fS,I] = sort(f0(:)); f0(I) = linspace(0,1,length(f0(:)));

% Laplacian
a = [2:n 1]; b = [n 1:n-1]; % per
a = [2:n n]; b = [1 1:n-1]; % sym
% grad
Grad = @(f)cat(3, f(a,:)-f, f(:,a)-f)/2;
Div = @(v)( v(b,:,1)-v(:,:,1) + v(:,b,2)-v(:,:,2) )/2;
Delta = @(f)-(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

switch eqname
    case 'heat'
        niter = 600;
        tau = .5;
    case 'tv'
        niter = 5000*2*4;
        epsilon = .001/4;
        tau = .5*epsilon;
end



q = 70;
k = 1;
ndisp = round(linspace(1,niter,q));
t = linspace(0,1,n);
 
f  = f0;
for i=1:niter    
    if i==ndisp(k)
        u = (i-1)/(niter-1);
        if 0
        % display 
        g = f/max(abs(f(:)));
        %
        r = 12; % #levellines
        clf; hold on;
        imagesc(t,t,g);
        contour(t,t,g,linspace(0,1,r), 'k');
        M = linspace(0,1,r-1)';
        colormap(M*[u 0 1-u] + (1-M)*[1 1 1]);
        caxis([0 1]);
        else
            clf; imagesc(t,t,f); colormap gray(256);
        end
        axis image; axis off; axis xy;
        drawnow;
        % saveas(gcf, [rep 'evol-' znum2str(k,3) '.png'], 'png');
        imwrite(rescale(f), [rep 'evol-' znum2str(k,3) '.png']);
        k = k+1;
    end
    switch eqname
        case 'heat'
            f = f - tau * Div( Grad(f) );
        case 'tv'
            G = Grad(f);
            G = G ./ repmat( sqrt(epsilon^2 + sum(G.^2,3)), [1 1 2] );
            f = f - tau * Div( G );            
    end
end