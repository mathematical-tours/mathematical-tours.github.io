% test for different cuves 

addpath('../toolbox/');
rep = MkResRep();

n = 2048;

p = 9; % #curvs
rho = linspace(0,1,p);

xmax = 5;
t = linspace(-pi,pi,n);
u = t+sin(t)+pi;
v = (1+cos(t));
x = linspace(0,xmax,n)';
y0 = interp1(u,v,x);


curve = @(rho)(2*rho-1)*y0(end)/x(end)*x + (1-(2*rho-1))*y0;
%
curve = @(rho)1*(2*rho-1)*sin(2*pi*x/xmax) + y0;
%
curve = @(rho).8*(2*rho-1)*sin(4*pi*x/xmax) + y0;


% Brachistochrone equation

T = [];
Y = []; 
Yd = [];
for j=1:p
    Y(:,j) = curve(rho(j));
    Yd(:,j) = n/xmax * diff(Y(:,j)); 
    % t as a function of x
    U(:,j) = sqrt((1+Yd(:,j).^2)./Y(2:end,j));
    T(:,j) = xmax/n * [0; cumsum(U(:,j))];
end

clf;
plot(x,1-Y);

L = max(T(:));
q = 80; % #times
t = linspace(0,L,q);

Xt = [];Yt = [];
for j=1:p
    Xt(:,j) = interp1(T(:,j), x, t);
end
Xt(isnan(Xt))=xmax;
for j=1:p
    Yt(:,j) = interp1(x, Y(:,j), Xt(:,j));
end

% animate
for it=1:q
    clf; hold on;
    for j=1:p
        lw = 2;
        ms = 20;
        if j==(p+1)/2
            lw = 4;
        end
        s = (j-1)/(p-1);
        plot(x,1-Y(:,j), 'Color', [s 0 1-s], 'LineWidth', lw);
        plot(Xt(it,j),1-Yt(it,j), '.', 'MarkerSize', ms, 'Color', [s 0 1-s]);
    end
    axis off;
    axis equal;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png'] );
end
