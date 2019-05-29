%%
% Non linear heat equation.


H = 2;
H = 1;
H = 3;
H = 4;


addpath('../toolbox/');
rep = MkResRep(['1d-' num2str(H)]);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% size
n = 256*2;

% Laplacian
a = [2:n n]; b = [1 1:n-1];
a = [2:n 1]; b = [n 1:n-1]; % per
Delta = @(f)-(2*f - f(a) - f(b) )/4;

% initialization
t = linspace(0,1,n)';

k = 2;
rand('state', 123);
Z = rand(k,1);
w = 1 + .2*rand(k,1);

Z = [.2 .6 .8];
w = [.5 1 .8];


s = .02;
f0 = zeros(n,1);
for k=1:length(Z)
    z = Z(k);
    f0 = f0 + w(k) * exp(-( (t-z).^2 )/(2*s^2));
end



switch H
    case 4
        Tmax = 20000*3;
        tau = .9;
    case 3
        Tmax = 25000*2;
        tau = .6;
    case 2
        Tmax = 10000;
        tau = .9;
    case 1
        Tmax = 3000;
        tau = .8;
end
niter = round(Tmax/tau);


% svg run
q = 50; 
ndisp = round(linspace(1,niter,q));
k = 1; 
%
f  = f0;
f1 = f0;
F = zeros(n,q);
for i=1:niter  
    if ndisp(k)==i
        F(:,k) = f/sum(f);
        k = k+1;
    end
    f = f + tau * Delta(f.^H);
end
axis tight;


% AutoCrop(rep, ['porous-']);
% convert heat-*.png heat.gif

% re-render
for i=1:q
    clf; hold on;
    for j=1:i
        u = (j-1)/(q-1);
        plot(t, F(:,j), 'Color', [1 1 1]*.7 + .3*[u 0 1-u], 'LineWidth', 1);
    end
    u = (i-1)/(q-1);
    plot(t, F(:,i), 'Color', [u 0 1-u], 'LineWidth', 2);
    axis([0 1 0 max(F(:))*1.01]);
    % axis([0 1 -.02 1.02]);
    box on;
    set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/2);
    drawnow;     
    saveas(gcf, [rep  'anim-' znum2str(i,2) '.png'], 'png');
end

%  AutoCrop(rep, 'anim')