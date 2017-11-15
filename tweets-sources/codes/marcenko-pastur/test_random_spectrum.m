%%
% test for various random matrix distributions

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

rep = ['results/randmat/'];
[~,~] = mkdir(rep);


ms = 20;
n = 1000;

t = linspace(0,1,1000);


SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

%% 
% randn

A = randn(n)/sqrt(n);
clf; hold on;
plot( eig(A), '.', 'MarkerSize', ms  );
plot( cos(2*pi*t), sin(2*pi*t), 'r', 'LineWidth', 2 )
axis tight;
axis equal; box on;
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'circlelaw.eps'], 'epsc');

%% 
% quarter circle

A = randn(n)/sqrt(n);
A = (A+A')/sqrt(2);
A = A-diag(diag(A))+diag(randn(n,1)/sqrt(n));
u = eig(A);


q = 60;
L = 2.3;
r = linspace(-L,L,q);
h = hist(u,r); h = q*h/(2*L)/sum(h);

Q = 1000;
T = linspace(-L,L,Q);
%
QC = sqrt( 4-T.^2  );
QC( T<-2 | T>2 ) = 0;
QC = Q*QC/(2*L)/sum(QC);


%
clf; hold on;
bar(r,h);
plot(T,QC, 'r', 'LineWidth', 2);
axis tight; box on;
set(gca, 'FontSize', 20);
SetAR(1/2);
saveas(gcf, [rep 'quartercircle.eps'], 'epsc');


%% 
% marcenko_pastur

beta = 1/10;
beta = .9;
p = round(n/beta);
B = randn(p,n)/sqrt(p);
A = B'*B;
u = eig(A);

Q = 1000;
T = linspace(1e-5,2,Q);
%
a = (1-sqrt(beta))^2;
b = (1+sqrt(beta))^2;
MP = sqrt((b-T).*(T-a))./( 2*pi*beta*T );
MP( T<a | T>b ) = 0;
MP = Q*MP/sum(MP);

q = 60;
r = linspace(0,2,q);
h = hist(u,r); h = q*h/sum(h);
%
clf; hold on;
bar(r,h);
plot(T,MP, 'r', 'LineWidth', 2);
axis tight; box on;
set(gca, 'FontSize', 20);
SetAR(1/2);
saveas(gcf, [rep 'marcenko_pastur.eps'], 'epsc');
