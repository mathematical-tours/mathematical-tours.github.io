%%
% Test Monte Carlo

d = 3;  % #dimension

addpath('../toolbox/');
rep = MkResRep(num2str(d));

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


%%
% Just a display.

nmax = 1000;
q = 50; % for anim

ndisp = round(linspace(1,nmax,q));

X = rand(d,nmax);

r = 1/2;
I = sum((X-1/2).^2)<=r^2;
J = sum((X-1/2).^2)>r^2;


for it=1:q
    s = (it-1)/(q-1);
    n = ndisp(it);
    
    clf;  hold on;
    switch d
        case 3
            k = 100;
            [x,y,z] = sphere(k);
            surf(1/2+x/2,1/2+y/2,1/2+z/2, ones(k+1)); alpha(.2);
            shading interp;
            plot3(X(1,I(1:n)),X(2,I(1:n)),X(3,I(1:n)), '.', 'Color', [s 0 1-s], 'MarkerSize', 25);
            plot3(X(1,J(1:n)),X(2,J(1:n)),X(3,J(1:n)), '.',  'Color', .8+.2*[s 0 1-s], 'MarkerSize', 25);
            view(3);
            axis equal; axis on; box on;
            set(gca, 'Xtick', [], 'Ytick', [], 'Ztick', []);
            camlight();
        case 2
            t = linspace(0,1,1000)*2*pi;
            plot(1/2+cos(t)/2, 1/2+sin(t)/2, 'g', 'LineWidth', 2);
            plot(X(1,I(1:n)),X(2,I(1:n)), '.', 'Color', [s 0 1-s], 'MarkerSize', 25);
            plot(X(1,J(1:n)),X(2,J(1:n)), '.',  'Color', .8+.2*[s 0 1-s], 'MarkerSize', 25);
            axis equal; axis([0 1 0 1]); axis on; box on; set(gca, 'Xtick', [], 'Ytick', [], 'Ztick', []);
    end
    drawnow;
    saveas(gcf, [rep 'points-' znum2str(it,2) '.png']);    
end
%

% target vol
switch d
    case 1
        V = 2*r;
    case 2
        V = pi*r^2;
    case 3
        V = (4/3)*pi*r^3;
end


% several simulations
nrep = 100;
VM = [];
for k=1:nrep
    X = rand(d,n);
    I = sum((X-1/2).^2)<=r^2;
    VM(:,end+1) = cumsum(I(:)) ./ (1:n)';
end


% plot evolution
K = 10;
C = distinguishable_colors(K);

for it=1:q
    s = (it-1)/(q-1);
    n = ndisp(it);
    
    clf; hold on;
    for k=2:K
        plot(VM(:,k), 'Color', [1 1 1]*.7, 'LineWidth', 1);
    end
    plot(VM(:,k), 'Color', [0 0 0], 'LineWidth', 3);
    plot(n,VM(n,k), '.', 'Color', [s 0 1-s], 'MarkerSize', 45);
    plot([1 nmax], [V V], 'k--', 'LineWidth', 2);
    plot(1:nmax, V+1./sqrt(1:nmax), 'k--', 'LineWidth', 2);
    plot(1:nmax, V-1./sqrt(1:nmax), 'k--', 'LineWidth', 2);
    mybox([1 nmax .7*V 1.3*V], [s 0 1-s], 3);
    axis([1 nmax .7*V 1.3*V]);
    set(gca, 'XTick', [], 'YTick', []);
    box on;
    SetAR(1/2);
    set(gca, 'FontSize', 20);
    drawnow;
    saveas(gcf, [rep 'integral-' znum2str(it,2) '.png']);   
    
end

return;

E = mean(abs(VM-V), 2);
clf; 
plot(log10(1:n),log10(E), 'LineWidth', 2);
SetAR(1/2);
axis tight;

% AutoCrop(rep, 'points');
% AutoCrop(rep, 'integral');