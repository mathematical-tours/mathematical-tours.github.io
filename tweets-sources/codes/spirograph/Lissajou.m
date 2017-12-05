
rep = 'results/lissajou/';
[~,~] = mkdir(rep); 

rho = .05;
t = 2*pi*linspace(0,1/rho,2048);

q = 1+rho;


Q = 100;

clf; hold on;
for i=1:Q
    c = (i-1)/Q;
    t = linspace( (i-1)/Q, i/Q,100 )*1/rho*2*pi;
plot(sin(t),cos((1+rho)*t), 'LineWidth', 2, 'Color', [c 0 1-c]);
axis([-1 1 -1 1]);
axis off;
drawnow;
    saveas(gcf, [rep 'liss-' num2str(i) '.png'], 'png');
end