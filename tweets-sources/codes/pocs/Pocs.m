%%
% Iterative projection on convex and non-convex sets.

addpath('./toolbox');

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
name = 'sines';

% Select which algorithm POCS or DR
algotype = 'proj';
algotype = 'dr';

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
    case 'sines'
        t = linspace(-1,1,5000)';
        c = 5*t + .5*1i*cos(t*5*pi);
        c1 = c*exp(2i*pi*.12);
        AX = [-3 3 -3 3];
        Proj = { @(x)ProjCurve(x,c),@(x)ProjCurve(x,c1) };
        Draw = { @()c,@(x)c1 };
end


niter = 20;
ntrials = 4;
C = distinguishable_colors(ntrials+4); C = C(5:end,:);

% reflexive projection
RProj = {@(x)2*Proj{1}(x)-x, @(x)2*Proj{2}(x)-x};

if not(exist('a'))
    a = {};
end

ms = 15; % marker size
clf; hold on;
plot( Draw{1}(), 'b', 'LineWidth', 2);
plot( Draw{2}(), 'r', 'LineWidth', 2);
axis(AX);
for it=1:ntrials
    x = randn + 1i*randn;
    % for DR
    mu = 1;
    if length(a)<it
        [a{it},b{it}] = ginput(1); 
    end
    x = a{it}+1i*b{it};
    tx = x;
    for i=1:niter
        switch algotype
            case 'proj'
                x(end+1) = Proj{1}(x(end));
                x(end+1) = Proj{2}(x(end));
            case 'dr'
                tx = (1-mu/2)*tx + mu/2*RProj{2}( RProj{1}(tx) );
                x(end+1) = RProj{1}(tx);
        end
    end
    plot(x, '.-', 'LineWidth', 1, 'MarkerSize', ms, 'color', C(it,:));
end
axis equal;
axis(AX);
axis off;
saveas(gcf, [rep name '-' algotype '.eps'], 'epsc');
