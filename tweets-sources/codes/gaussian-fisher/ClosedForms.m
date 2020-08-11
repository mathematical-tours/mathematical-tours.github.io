%%
% Closed form for divergences between Gaussians

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep 'a-' name '-' znum2str(it,3) '.png']);


%  'wasserstein'
W = @(m1,s1,m2,s2)sqrt( (m1-m2).^2 + (s1-s2).^2 );
% 'hellinger'
A = @(m1,s1,m2,s2)sqrt(2*s1.*s2./(s1.^2+s2.^2)) .* ...
    exp( -(m1-m2).^2./( 2*(s1.^2+s2.^2) ) );
H = @(m1,s1,m2,s2)sqrt(2-2*A(m1,s1,m2,s2));
% 'kl'
KL = @(m1,s1,m2,s2)log(s2./s1) + (s1.^2 + (m1-m2).^2)./(2*s2.^2) - 1/2;
% Fisher
F = @(m1,s1,m2,s2)( (m1-m2).^2 + 2*(s1-s2).^2 ).* ...
        ( (m1-m2).^2 + 2*(s1+s2).^2 );
F = @(m1,s1,m2,s2)sqrt( F(m1,s1,m2,s2) );
Fisher = @(m1,s1,m2,s2)sqrt(2)*log( ...
        ( F(m1,s1,m2,s2) + (m1-m2).^2 + 2*(s1.^2+s2.^2) ) ./ (4*s1.*s2) ...
    );

m1 = 0; m2 = 1; 
s1 = 1; s2 = 1.2;
% Fisher(m1,s1,m2,s2) - Fisher(m2,s2,m1,s1)
% H(m1,s1,m2,s2) - H(m2,s2,m1,s1)

mN = 300; sN = 200;
mN = 1000; sN = 600;
m = linspace(-1,1,mN); s = linspace(1e-3,1.5,sN);
[S,M] = meshgrid(s,m);

if not(exist('name'))
name = 'wasserstein';
name = 'klsym';
name = 'kl';
name = 'burg';
name = 'hellinger';
name = 'fisher';
end

switch name
    case 'wasserstein'
        f = W; vmax = Inf;
    case 'hellinger'
        f = H; vmax = Inf;
    case 'kl'
        f = @(m1,s1,m2,s2)sqrt(KL(m1,s1,m2,s2)); vmax = 10;
    case 'burg'
        f = @(m1,s1,m2,s2)sqrt(KL(m2,s2,m1,s1)); vmax = 8;
    case 'fisher'
        f = @(m1,s1,m2,s2)sqrt(Fisher(m1,s1,m2,s2)); vmax = 5;
    case 'klsym'
        f = @(m1,s1,m2,s2)sqrt(KL(m2,s2,m1,s1) + KL(m1,s1,m2,s2)); vmax = 10;
end

sa = .75; 
q = 50;
for it=1:q
    t = (it-1)/q;
    m1 = cos(2*pi*t)*.7;
    s1 = sa + sin(2*pi*t)*sa*.8;
    
    
    U = f(m1,s1,M,S);
    U = min(U,vmax);
    % U = min(U,max(U(:,50)));
    U = U/max(U(:));
    
    r = 18; % #levellines
    clf; hold on;
    imagesc(m,s,U');
    contour(m,s,U',linspace(0,1,r), 'k');
    % contour(m,s,f(m1,s1,M,S)',f(m1,s1,m1*ones(1,r),linspace(1e-2,s1,r)), 'k');
    plot(m1,s1, 'k.', 'MarkerSize', 20);
    colormap(parula(r-1));
    % colormap(parula(256));
    caxis([0 1]);
    axis equal; axis off;
    drawnow;
    mysaveas(name,it);
end


