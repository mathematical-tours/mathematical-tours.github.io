%%
% simple 2D example of Gibbs sampler

addpath('../toolbox/');
rep = MkResRep();

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
caxis([0 1]);

% temperature
eps = .9;
% movement noise
sigma = .01/30;
niter = 100000;


% temperature
eps = .2;
% movement noise
sigma = .01/6;
niter = 40000;

% temperature
eps = .1/2;
% movement noise
sigma = .01/5;
niter = 4000;




% movement noise
eps = .1/10;
sigma = .01/2;
niter = 500;

% movement noise
eps = .1/2;
sigma = .01/2;
niter = 900;

% display the gibbs kernel
G = exp(-f(X+1i*Y)/eps);
r = 14;
clf; hold on;
imagesc(t,t,G');
contour(t,t,G',linspace(0,1,r), 'k');
colormap(parula(r-1));
caxis([0 1]);

% monte-carlo Gibb sampling
p = 40; % #particles
m = 60; % memory for display
% initial seed
z = 1/2+1i/2 + (randn(p,1)+1i*randn(p,1))*.01;


q = 90; 
ndisp = round(linspace(1,niter,q));
kdisp = 1;

fast_mode=0;

for it=1:niter
    % old point
    z0 = z(:,end);
    % new proposals
    z1 = z0 + (randn(p,1)+1i*randn(p,1))*sigma;
    % accept/reject
    R = exp((f(z0)-f(z1))/eps);
    U = rand(p,1)>R;
    z(:,end+1) = z1.*(1-U) + z0.*U;
   
    % erase old points
    z(:,1:end-m) = [];
    
    if ndisp(kdisp)==it
    clf; hold on;
    imagesc(t,t,F');
    contour(t,t,F',linspace(0,1,r), 'k');
    colormap(parula(r-1));
    caxis([0 1]);
    % display trajectories
    if fast_mode==1
        plot(z(:,end), 'k.', 'LineWidth', 2, 'MarkerSize', 25);
    else
        for j=1:size(z,2)-1
            s = (j-1)/(size(z,2)-1);
            set(gca,'ColorOrderIndex',1)
            h = plot( real(z(:,j:j+1))', imag(z(:,j:j+1))', 'LineWidth', 3); % , 'Color', [1 0 0] );
            for k=1:p % transparency
                h(k).Color(4) = s;
            end
        end
        % set(gca,'ColorOrderIndex',1)
        % plot(z(:,end), '.', 'MarkerSize', 25);
        %scatter(real(z(:,end)), imag(z(:,end)));
    end
    axis([0 1 0 1]); axis off;
    drawnow;
    if fast_mode==0
        saveas(gcf, [rep 'anim-' znum2str(kdisp,2) '.png']);
    end
    kdisp = kdisp+1;
    end
end

% AutoCrop(rep, 'anim');