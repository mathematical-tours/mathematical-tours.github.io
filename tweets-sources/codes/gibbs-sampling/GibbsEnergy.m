%% 
% Just display Gibbs energy and sampling


addpath('../toolbox/');
rep = MkResRep('gibbs');

n = 200; % for display
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

% define a function
q = 1.3;
gauss = @(z,c,s)exp(- abs(z-c).^q/(2*s^q));
c = [.3 + .7i, .6 + .25i, .65 + .75i]; 
s = [.1, .13, .08];
f = @(z)1 - .8*gauss(z,c(1),s(1)) - .9*gauss(z,c(2),s(2)) - gauss(z,c(3),s(3)) ;


% display the function
F = f(X+1i*Y);
r = 14;
clf; hold on;
imagesc(t,t,F');
contour(t,t,F',linspace(0,1,r), 'k');
colormap(parula(r-1));
caxis([0 1]); axis off;
saveas(gcf, [rep 'func-' znum2str(k,2) '.png']);

m = 100; % # samples

% temperature
q = 10;
eps_list = linspace(.8,0.05,q);

k = 1;
for i=1:q
    eps = eps_list(i);
    P = exp(-F/eps); P = P/sum(P(:));
    
    I = rand_discr(P(:),m);
    z = X(I)+1i*Y(I); 

    % colormap interpolating between two colors (here a color and white)
    r = 14;
    s = (i-1)/(q-1);
    h = linspace(0,1,r-1)';
    CM = h*[s 0 1-s] + (1-h)*[1 1 1];

    P = P/max(P(:));
    clf; hold on;
    imagesc(t,t,P');
    contour(t,t,P',linspace(0,1,r), 'k');
    plot(z, 'k.', 'MarkerSize', 25);
    colormap(CM);
    caxis([0 1]); axis off;
    drawnow;
    for j=1:5
        saveas(gcf, [rep 'anim-' znum2str(k,2) '.png']);
        k = k+1;
    end
end

% AutoCrop(rep, 'anim');