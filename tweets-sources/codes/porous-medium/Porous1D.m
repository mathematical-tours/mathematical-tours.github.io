%%
% Non linear heat equation.


addpath('../toolbox/');
rep = MkResRep('1d');
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



H = 3;
H = 1;
switch H
    case 4
        Tmax = 200000;
        tau = .9;
    case 3
        Tmax = 25000;
        tau = .6;
    case 2
        Tmax = 10000;
        tau = .9;
    case 1
        Tmax = 4000;
        tau = .8;
end
niter = round(Tmax/tau);



% svg run
ndisp = max(1,ceil(niter/50)); 
k = 0; 
%
f  = f0;
f1 = f0;
for i=1:niter    
    if mod(i,ndisp)==1
        k = k+1;
        u = (i-1)/(niter-1);
        % display 
        g = f / max(abs(f(:)));
        clf;
        plot(t, g, 'Color', [u 0 1-u], 'LineWidth', 2);
        axis([0 1 -.02 1.02]);
        box on;
        set(gca, 'XTick', [], 'YTick', []);
        SetAR(1/2);
        drawnow;
        %
        saveas(gcf, [rep  'porous-' num2str(H) '-' znum2str(k,2) '.png'], 'png');
    end
    f = f + tau * Delta(f.^H);
end
axis tight;


% AutoCrop(rep, ['porous-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif