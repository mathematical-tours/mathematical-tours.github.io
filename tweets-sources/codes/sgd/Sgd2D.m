%%
% Test for SGD in 2D.

name = 'aniso';
name = 'iso';

addpath('../toolbox/');
rep = MkResRep(name);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

switch name
    case 'iso'
        eta = 1;
        vwarm = 10;
        tau0 = 1;
    case 'aniso'
        eta = .4;
        vwarm = 30;
        tau0 = 3;
end

r=20; % #particles
crandn = @(k)randn(k,1)+1i*randn(k,1);
x = [crandn(r/2)*.3-1/2; crandn(r/2)*.3+1/2];

x = cos( (0:r-1)'/r*2*pi ) + 1i * sin( eps + (0:r-1)'/r*2*pi );

% int E( <x,theta>^2 )

niter = 1000;
E = [];
for i=1:niter-1
    u = randn(r,1) + 1i*randn(r,1)*eta;
    E(:,i) = abs(x(:,end)).^2;
    tau = tau0/(vwarm+i);
    x(:,end+1) = x(:,end) - tau * real(conj(u).*x(:,end)) .* u;
end

% #frames
q = 80; 
C = distinguishable_colors(r);

% background image
N = 201;
ts = linspace(-1,1,N);
[Y,X] = meshgrid(ts,ts);
F = sqrt(X.^2+(eta*Y).^2);
rl = 15; % #levellines

for i=1:q
    t = (i-1)/q;
    j = 1+ceil( t^3 * niter );
    j = max(j,i);
    clf; hold on;
    % bagkground
    imagesc(ts,ts,F');
    contour(ts,ts,F',linspace(0,max(F(:)),rl), 'k');
    colormap(parula(rl-1));
    caxis([0 max(F(:))]);
    % curves
    for m=1:r
        plot(transpose(x(m,[j j])), '.', 'MarkerSize', 25, 'color', C(m,:));
        if j>1
            plot(transpose(x(m,1:j)), '-', 'LineWidth', 2, 'color', C(m,:));
        end
    end
    axis off;
    axis equal; axis([-1 1 -1 1]*1);
    drawnow;
    saveas(gcf, [rep 'sgd-' znum2str(i,2) '.png']);
end


% AutoCrop(rep, ['sgd-']); 

