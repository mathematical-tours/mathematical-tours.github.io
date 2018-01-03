%%
% Iterative projection on convex and non-convex sets.


rep = '../results/pocs/';
[~,~] = mkdir(rep);

t = linspace(0,1,200);
%
ProjCirc = @(x,c,r)c + r * (x-c)/abs(x-c);
DrawCirc = @(c,r)c + exp(2i*pi*t)*r;
% <x,u>=0
ProjLine = @(x,u)x - real(x*conj(u))*u/abs(u)^2;
DrawLine = @(u)20*(t-.5)*( imag(u)-1i*real(u) );

name = 'linecircle';
name = 'twolines';
name = 'twocircles';

AX = [-2 2 -2 2];
switch name
    case 'linecircle'
        u = exp(1i*pi/3);
        Proj = { @(x)ProjLine(x,u),@(x)ProjCirc(x,0.8,1) };
        Draw = { @()DrawLine(u),@(x)DrawCirc(0.8,1) };
        AX = [-1 2.5 -2 1.5];
    case 'twolines'
        u = exp(1i*pi/3);
        v = exp(-1i*pi/3);
        Proj = { @(x)ProjLine(x,u),@(x)ProjLine(x,v) };
        Draw = { @()DrawLine(u),@(x)DrawLine(v) };
    case 'twocircles'
        Proj = { @(x)ProjCirc(x,.5,1),@(x)ProjCirc(x,-.5,1) };
        Draw = { @()DrawCirc(.5,1),@(x)DrawCirc(-.5,1) };
end


niter = 20;
ntrials = 4;
C = distinguishable_colors(ntrials+4); C = C(5:end,:);

clf; hold on;
plot( Draw{1}(), 'b', 'LineWidth', 2);
plot( Draw{2}(), 'r', 'LineWidth', 2);
axis(AX);
for it=1:ntrials
    x = randn + 1i*randn;
    [a,b] = ginput(1); x = a+1i*b;
    for i=1:niter
        x(end+1) = Proj{1}(x(end));
        x(end+1) = Proj{2}(x(end));
    end
    plot(x, '.-', 'LineWidth', 1, 'MarkerSize', 25, 'color', C(it,:));
end
axis equal;
axis(AX);
axis off;
saveas(gcf, [rep name '.eps'], 'epsc');
