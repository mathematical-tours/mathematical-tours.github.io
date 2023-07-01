%%
% Display progressive power series.

addpath('../toolbox/');
rep = MkResRep();

dmax = 30;

% series
a = (-1).^(0:dmax-1);
f0 = @(t)1./(1+t);

%
a = (-1).^(0:dmax-1);
f0 = @(t)1./(1+t);
%
a = (+1).^(0:dmax-1);
f0 = @(t)1./(1-t);
xrange = [-1.5 1.1];
yrange = [-.5 4];



%
dmax = 15;
a = 1 ./ factorial(0:dmax-1);
f0 = @(t)exp(t);
xrange = [-5 2];
yrange = [-2 4];




%
dmax = 15;
a = (1i).^(0:dmax-1) ./ factorial(0:dmax-1);
f0 = @(t)exp(1i*t);
xrange = [-10 10];
yrange = [-2 2];

%
a = zeros(dmax,1);
a(2:2:end) = (-1).^(0:dmax/2-1) ./ factorial( 1:2:dmax-1 ); 
f0 = @(t)sin(t);
xrange = 12*[-1 1];
yrange = [-1.8 1.8];


%
dmax = 15*4;
a = zeros(dmax-1,1); 
a(1:2:end) = (-1).^(0:dmax/2-1);
f0 = @(t)1./(1+t.^2);
xrange = [-2.5 2.5];
yrange = [-.5 2.5];


% grid
n = 501;
N = 2048;
t = linspace(xrange(1),xrange(2),N)';  %for 1D
T = linspace(xrange(1),xrange(2),n)';
[Y,X] = meshgrid(T,T);
Z = X+1i*Y;



F = zeros(n,n);
f = zeros(N,1);
S = 10; % sub-animations
S = 5;
it = 0;
Ideg = find(a~=0);
for d=Ideg(:)'
    F1 = F + a(d) .* Z.^(d-1);
    f1 = f + a(d) .* t.^(d-1);
    for s=1:S
        it = it+1;
        u = (it-1)/(S*length(Ideg)-1);
        Fd = (1-s/S)*F + s/S * F1;
        fd = (1-s/S)*f + s/S * f1;
        %
        clf; hold on;
        if 1
            if isreal(f0(t))
                plot(t,f0(t), 'k--', 'LineWidth', 1);
                plot(t,fd, 'color',[u 0 1-u], 'LineWidth', 2);
                axis([xrange(1) xrange(2) yrange(1) yrange(2)]);
                set(gca, 'PlotBoxAspectRatio', [1 2/3 1]); 
            else
                plot(f0(t), 'k--', 'LineWidth', 1);
                plot(fd, 'color',[u 0 1-u], 'LineWidth', 2);
                axis square;
                axis([yrange(1) yrange(2) yrange(1) yrange(2)]);
            end
            %
        else
            PolyDisp(Fd,T);
        end
        box on;           
        set(gca, 'FontSize', 20, 'XTick', [], 'YTick', []);
        drawnow;
        saveas(gcf, [rep  'anim-' znum2str(it,3) '.png'] );
    end
    F = F1; f = f1;
end