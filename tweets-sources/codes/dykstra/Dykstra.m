%%
% Iterative projection vs Dykstra methods.

addpath('../toolbox');

rep = '../results/dykstra/';
[~,~] = mkdir(rep);

t = linspace(0,1,200);
% |x-c|=r
ProjCirc = @(x,c,r)c + r * (x-c)/abs(x-c);
DrawCirc = @(c,r)c + exp(2i*pi*t)*r;
% |x-c|<=r
DrawDisk = @(c,r)DrawCirc(c,r);
% <x,u>=0
ProjLine = @(x,u)x - real(x*conj(u))*u/abs(u)^2;
DrawLine = @(u)20*(t-.5)*( imag(u)-1i*real(u) );
% <x,u><=0
DrawHalf = @(u)DrawLine(u);


name = 'linecircle';
name = 'twocircles';
name = 'sines';
name = 'halfdisk';
name = 'linedisk';
name = 'twohalf';
name = 'halfline';

% Select which algorithm POCS or DR
algotype = 'dr';
algotype = 'dykstra';
algotype = 'proj';

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
    case 'twohalf'
        u = exp(1i*pi/3);
        v = -exp(-1i*pi/3);
        Proj = { @(x)ProjHalf(x,u),@(x)ProjHalf(x,v) };
        Draw = { @()DrawLine(u),@(x)DrawLine(v) };
    case 'halfline'
        u = exp(1i*pi/3);
        v = -exp(-1i*pi/3);
        Proj = { @(x)ProjHalf(x,u),@(x)ProjLine(x,v) };
        Draw = { @()DrawHalf(u),@(x)DrawLine(v) };
    case 'linedisk'
        u = exp(1i*pi/3);
        c = 0;
        Proj = { @(x)ProjDisk(x,c,1), @(x)ProjLine(x,u) };
        Draw = { @(x)DrawDisk(c,1), @()DrawLine(u) };
        AX = [-1 2.5 -2 1.5];
        AX = [-1 1 -1 1]*3;
    case 'halfdisk'
        u = exp(1i*pi/3);
        c = 0;
        Proj = { @(x)ProjDisk(x,c,1), @(x)ProjHalf(x,u) };
        Draw = { @(x)DrawDisk(c,1), @()DrawHalf(u) };
        AX = [-1 2.5 -2 1.5];
        AX = [-1 1 -1 1]*3;
    case 'twocircles'
        Proj = { @(x)ProjCirc(x,.5,1),@(x)ProjCirc(x,-.5,1) };
        Draw = { @()DrawCirc(.5,1),@(x)DrawCirc(-.5,1) };
    case 'sines'
        t = linspace(-1,1,5000)';
        c = 5*t + .5*1i*cos(t*5*pi);
        c1 = c*exp(2i*pi*.12);
        AX = [-3 3 -3 3];
        Proj = { @(x)ProjCurve(x,c),@(x)ProjCurve(x,c1) };
        Draw = { @()c,@(x)c1 };
end

% flip
if 0
Proj = {Proj{2} Proj{1}};
Draw = {Draw{2} Draw{1}};
end

niter = 20;
ntrials = 6;
C = distinguishable_colors(ntrials+4); C = C(5:end,:);

% reflexive projection
RProj = {@(x)2*Proj{1}(x)-x, @(x)2*Proj{2}(x)-x};

if not(exist('a'))
    a = {};
end

% for DR
mu = 1;
    
ms = 20; % marker size
clf; hold on;
plot( Draw{1}(), 'b', 'LineWidth', 2);
plot( Draw{2}(), 'r', 'LineWidth', 2);
axis(AX);
for it=1:ntrials
    if length(a)<it
        axis equal;
        [a{it},b{it}] = ginput(1); 
    end
    x = a{it}+1i*b{it};
    % for dr
    tx = x; 
    % for dkstra
    p=0; q=0;
    for i=1:niter
        switch algotype
            case 'proj'
                x(end+1) = Proj{1}(x(end));
                x(end+1) = Proj{2}(x(end));
            case 'dr'
                tx = (1-mu/2)*tx + mu/2*RProj{2}( RProj{1}(tx) );
                x(end+1) = RProj{1}(tx);
            case 'dykstra'
%                 y = Proj{1}(x(end) + p);
%                 p = x(end) + p - y;
%                 x(end+1) = Proj{2}(y+q);
%                 q = y+q-x(end);
                y = Proj{1}(x(end) + p);
                p = x(end) + p - y;
                x1 = Proj{2}(y+q);
                q = y+q-x1;
                x(end+1) = y; x(end+1) = x1;
        end
    end
    plot(x, '.-', 'LineWidth', 1, 'MarkerSize', ms, 'color', C(it,:));
    axis equal; axis(AX);
end
axis equal; axis(AX);
axis off;
saveas(gcf, [rep name '-' algotype '.eps'], 'epsc');
