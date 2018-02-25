%%
% Display of continuous and discrete sines

rep = '../results/fourier-atoms/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


%%
% 1D

% discrete
n = 32;
% continuous
N = 512;
% frequencies
fc = 4;
for f=1:fc    
    c = (f-1)/(fc-1);
    clf;
    plot(0:n, cos( 2*pi/n*f*(0:n) ), '.-', 'color', [c 0 1-c], 'LineWidth', 2, 'MarkerSize', 25);
    axis([0 n -1.05 1.05]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/4);
    saveas(gcf, [rep '1d-discr-' num2str(f) '.eps'], 'epsc');
    clf;
    plot(0:N, cos( 2*pi/N*f*(0:N) ), '-', 'color', [c 0 1-c], 'LineWidth', 2, 'MarkerSize', 25);
    axis([0 N -1.05 1.05]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/4);
    saveas(gcf, [rep '1d-cont-' num2str(f) '.eps'], 'epsc');
end

%% 
% 2D


% discrete
n = 16;
% continuous
N = 256;
% frequencies
Fx = 1:4;
Fy = round(1.5*Fx);

[y,x] = meshgrid(1:n,1:n);
[Y,X] = meshgrid(1:N,1:N);

S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));

for k=1:length(Fx) 
    c = (k-1)/(length(Fx)-1);
    c0 = [1 1 1]; c1 = [c 0 1-c];
    fx = Fx(k); fy = Fy(k);
    c = (f-1)/(fc-1);
    % discrete
    u = cos( 2*pi/n*( fx*x+fy*y ) );
    U = (1+Upsc(u,N/n))/2; % in [0,1]
    V = zeros(N,N,3);
    for r=1:3
        V(:,:,r) = (1-U)*c0(r) + U*c1(r);
    end
    imwrite(V, [rep '2d-discr-' num2str(k) '.png'], 'png');
    % continuous
    U = cos( 2*pi/N*( fx*X+fy*Y ) );
    U = (1+U)/2; % in [0,1]
    V = zeros(N,N,3);
    for r=1:3
        V(:,:,r) = (1-U)*c0(r) + U*c1(r);
    end
    imwrite(V, [rep '2d-cont-' num2str(k) '.png'], 'png');
end
