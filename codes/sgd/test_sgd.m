%%
% Test for SGD on a simple case of 2 functions

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

rep = ['./results/sgd/'];
[~,~] = mkdir(rep);

f = {@(x)1/2*(x-1).^2, @(x)1/2*(x+1).^2};
F = @(x)f{1}(x)+f{2}(x);
df = {@(x)(x-1), @(x)(x+1)};

q=10000; % #particles
x = rand(q,1)-1/2;

niter = 1000;
E = [];
for i=1:niter-1
    u = rand(q,1)>.5;
    E(:,i) = F(x(:,end));
    tau = 1/(10+i);
    x(:,end+1) = x(:,end) - tau * ( u.*df{1}(x(:,end)) + (1-u).*df{2}(x(:,end)) );
end


clf;
plot(log10(max(E'-F(0),1e-10)));
axis tight; box on; 

U = .5;
clf;
% plot(x(1:50,:)', 'r');
plot_semi_transparent(1:niter,x(1:200,:)','r',.1);
axis tight;
axis([1 niter -U U]);
SetAR(1/2); set(gca, 'FontSize', 20);  box on;
saveas(gcf, [rep 'sgd-trajectory.eps'], 'epsc');
saveas(gcf, [rep 'sgd-trajectory.png'], 'png');


% histogram 
r = 101;
t = linspace(-U,U,r);
H = hist(x,t);

clf; hold on;
imagesc(1:niter, t, -H); colormap gray(256);
axis tight; SetAR(1/2); set(gca, 'FontSize', 20); box on;
saveas(gcf, [rep 'sgd-histo.eps'], 'epsc');
saveas(gcf, [rep 'sgd-histo.png'], 'png');
