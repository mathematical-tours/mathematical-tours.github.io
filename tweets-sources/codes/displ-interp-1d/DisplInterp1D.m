%%
% Test for 1-D optimal transport.

N = 40;
vmax = 20;
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

rep = '../results/displ-interp-1d/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

% generate Gaussian mixtures points in 1-D

x = 14 + 2*randn(N,1);
x = round(clamp(x,0,vmax));
y = 6 + 3*randn(N,1);
y = round(clamp(y,0,vmax));

% compute "elevation" for disantengling
Mx = zeros(vmax+1,1); Ex = []; 
My = zeros(vmax+1,1); Ey = [];
for i=1:N
    Ex(i) = Mx(x(i)+1); Mx(x(i)+1) = Mx(x(i)+1)+1;
    Ey(i) = My(y(i)+1); My(y(i)+1) = My(y(i)+1)+1;
end
for i=1:N
    Ex(i) = Ex(i) - Mx(x(i)+1)/2;
    Ey(i) = Ey(i) - My(y(i)+1)/2;
end
Ex = Ex(:); Ey = Ey(:);

ms = 30;
r = .4;
% scatter plot
clf;
subplot(2,1,1);
plot(x, r*Ex, '.r', 'MarkerSize', ms); axis([0 vmax -.5 .5]); axis equal;
set(gca, 'YTick', []);
subplot(2,1,2);
plot(y, r*Ey, '.b', 'MarkerSize', ms); axis([0 vmax -.5 .5]); axis equal;
set(gca, 'YTick', []);


% histogram
subplot(2,1,1);
hx = hist(x,1:vmax);
bar(1:vmax,hx, 'r'); axis tight;
subplot(2,1,2);
hy = hist(y,1:vmax);
bar(1:vmax,hy, 'b'); axis tight;

%%
% animation on histograms

save_gif = 0;
K = 40;
% sort
[xs,Ix] = sort(x + 1e-3*Ex); ExI = Ex(Ix);
[ys,Iy] = sort(y + 1e-3*Ey); EyI = Ey(Iy);
%
for i=1:K
    %
    t = (i-1)/(K-1);
    z = (1-t)*xs+t*ys;
    Ez = (1-t)*ExI+t*EyI;
    hz = hist(z,1:vmax);
    % OT
    clf;
    plot(z, r*Ez, '.', 'color', [1-t,0,t], 'MarkerSize', ms); axis([0 vmax -2 2]); % axis equal;
    set(gca, 'YTick', [], 'XTick', [0:5:vmax], 'FontSize', 20);
    axis([0 vmax -3 3]);
    SetAR(1/3);
    saveas(gcf, [rep 'interp-ot-' num2str(i) '.png']);
    % hist OT
    clf;
    bar(1:vmax,hz, 'FaceColor', [1-t,0,t]); 
    M = max([hx(:);hy(:)]);
    axis([0,vmax,0,M*1.02]);
    set(gca, 'YTick', [0:2:M], 'XTick', [0:5:vmax], 'FontSize', 20);
    SetAR(1/2);
    saveas(gcf, [rep 'interp-hist-ot-' num2str(i) '.png']);
    % hist L2
    clf;
    bar(1:vmax,(1-t)*hx+t*hy, 'FaceColor', [1-t,0,t]); axis([0,vmax,0,8]);
    axis([0,vmax,0,M*1.02]);
    set(gca, 'YTick', [0:2:M], 'XTick', [0:5:vmax], 'FontSize', 20);
    SetAR(1/2);
    drawnow;
    saveas(gcf, [rep 'interp-hist-l2-' num2str(i) '.png']);
end

