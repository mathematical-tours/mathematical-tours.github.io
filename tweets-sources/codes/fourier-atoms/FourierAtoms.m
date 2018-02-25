%%
% Display of continuous and discrete sines

rep = '../results/fourier-atoms/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

%%
% 1D

% discrete
n = 32;
% continuous
N = 512;
% frequencies
fc = 8;
for f=1:fc    
    c = (f-1)/(fc-1);
    clf;
    plot(0:n, cos( 2*pi/n*(0:n) ), '.-', 'color', [c 0 1-c], 'LineWidth', 2, 'MarkerSize', 25);
    axis([0 n -1.05 1.05]); box on;
end