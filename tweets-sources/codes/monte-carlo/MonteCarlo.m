%%
%

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


%%
% Just a display.

n = 1000;
d = 3;  % #dimension
X = rand(d,n);

r = 1/2;
I = sum((X-1/2).^2)<=r^2;
J = sum((X-1/2).^2)>r^2;

clf;  hold on;
switch d
    case 3
        k = 100;
        [x,y,z] = sphere(k);
        surf(1/2+x/2,1/2+y/2,1/2+z/2, ones(k+1)); alpha(.2);
        shading interp;
        plot3(X(1,I),X(2,I),X(3,I), 'b.', 'MarkerSize', 25);
        plot3(X(1,J),X(2,J),X(3,J), 'r.', 'MarkerSize', 20);
        axis equal; axis on; box on; set(gca, 'Xtick', [], 'Ytick', [], 'Ztick', []);
        view(3);
        camlight();
    case 2
        t = linspace(0,1,1000)*2*pi;
        plot(1/2+cos(t)/2, 1/2+sin(t)/2, 'g', 'LineWidth', 2);
        plot(X(1,I),X(2,I), 'b.', 'MarkerSize', 25);
        plot(X(1,J),X(2,J), 'r.', 'MarkerSize', 20);
        axis equal; axis([0 1 0 1]); axis on; box on; set(gca, 'Xtick', [], 'Ytick', [], 'Ztick', []);
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
n = 5000;
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
clf; hold on;
for k=1:K
    plot(VM(:,k), 'Color', C(k,:), 'LineWidth', 1);
end
plot([1 n], [V V], 'r--', 'LineWidth', 2);
axis([1 n .8*V 1.2*V]); 
box on;
set(gca, 'FontSize', 20); 
SetAR(1/2);



E = mean(abs(VM-V), 2);
clf; 
plot(log10(1:n),log10(E), 'LineWidth', 2);
SetAR(1/2);
axis tight;