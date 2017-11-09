
addpath('toolbox_certif/');

name = 'gaussian1d';
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

rep = ['results/certif-convol/'];
[~,~] = mkdir(rep);

% number of points
N = 2;
N = 3;

% scaling of the problem
if not(exist('Delta'))
    Delta = 1;
end

% set to 1 to impose that the diagonal of the covariance is 1
normalize = 1;
[C,d,z0,xlim,ylim,opt.domain] = load_correlation(name,normalize);

%%
% Compute vanishing derivatives constraints

K = 2*N-1;
% list of polynomials describing the derivative constraints
clear P;
syms x;
for k=0:K
    P(k+1) = x^k;
end

switch N
    case 2
        x0 = [.3 .6]';
        s0 = [1 1]';
    case 3
        x0 = [.2 .5 .7]';
        s0 = [1 1 -1]';
end

x0 = 1/2 + (x0-1/2)*Delta;
X0 = [x0(:), x0(:)*0];
EtaV = compute_etav(C,X0,s0);
EtaF = compute_etaf(C,X0,s0);

Q = 800;
Xa = linspace(0,1,Q);
etaV = fast_formula(EtaV,Xa,Xa*0);
etaF = fast_formula(EtaF,Xa,Xa*0);


clf;
hold on;
plot(Xa,etaV, 'b', 'LineWidth', 2);
plot(Xa,etaF, 'r', 'LineWidth', 2);
plot([0 1], [1 1], 'k--');
plot([0 1],-[1 1], 'k--');
stem(x0,s0, 'k.', 'MarkerSize', 25);
axis tight; box on;
set(gca, 'FontSize', 15);
if min(s0)>0
    axis([0 1 -.1 1.1]);
else
    axis([0 1 -1.1 1.1]);
end
legend('\eta_V', '\eta_F');
SetAR(1/2);
saveas(gcf, [rep 'etaV-N' num2str(N) '-d' num2str(round(100*Delta)) '.eps'], 'epsc');
%% ZOOM
r = .02;
I = find( abs(Xa-1/2)<r );
axis([.5-r,.5+r,min(etaF(I)) 1.01*max(etaF(I))]);
saveas(gcf, [rep 'etaV-N' num2str(N) '-d' num2str(round(100*Delta)) '-Z.eps'], 'epsc');

  