%%
% test for Fourier transform, sampling, interpolation 

rep = 'results/';
[~,~] = mkdir(rep);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20, 'XTick', [], 'YTick', []);


name = 'hibiscus';
F = imread([name '.png']);

F = sum(F,3); F = (F-min(F(:)))/(max(F(:))-min(F(:)));
clf; imagesc(F); axis image; colormap gray;

imwrite(F, [rep 'img-original.png']);

% extract 1-D trace
f= F(:,end/2);
n = length(f);
clf; plot(f); axis tight;


% make it periodic 
t = linspace(0,1,n)';
f = f - f(1) + (f(1)-f(end))*t;
plot(f); axis tight;



% number of non-zero elements
Klist = 2*round( [1 .5 .2 .1]*n/2 ); Klist(end) = 4;
fh = fft(f);
for i=1:length(Klist)
    k = Klist(i);
    fh1 = zeros(n,1);
    fh1(1:k/2) = fh(1:k/2);
    fh1(end-k/2+1:end) = fh(end-k/2+1:end);
    f1 = real( ifft(fh1) );
    clf;
    plot(f1, 'LineWidth', 2); axis tight;
    SetAR(1/3);
    saveas(gcf, [rep 'signal-regul-' num2str(i) '.eps'], 'epsc');
end


% plot sinc
a = 6;
t = linspace(-a-.2,+a+.2,1024)+1e-15;
clf; hold on;
plot(t, sin(pi*t)./(pi*t), 'LineWidth', 2 );
plot(t, 0*t, 'k--',  'LineWidth', 1 );
plot(-a:a, zeros(2*a+1,1), 'r.', 'MarkerSize', 25);
axis tight; box on;
set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20);
saveas(gcf, [rep 'sinc.eps'], 'epsc');


% plot its fft
fh = fft(f);
plot( log(abs(fh)) ); axis tight;
plot( -n/2:n/2-1, log(abs(fftshift(fh))) ); axis tight;

% interpolate it though zero padding, i.e. sinc
q = 12; % interpolation factor
fh1 = zeros(q*n,1); 
fh1(1:n/2) = fh(1:n/2);
fh1(end-n/2+2:end) = fh(n/2+2:end);
fh1(n/2+1) = fh(end/2+1)/2; % fh(end/2+1) is the pivot point, needs to split it in 2
fh1(end-n/2+1) = fh(end/2+1)/2;
f1 = q * real( ifft(fh1) ); % must be real anyway


clf; hold on;
plot(1:q:n*q,f, 'b.', 'MarkerSize', 20);
plot(1:n*q, f1);
axis tight; box on;


% plot a subpart
sel = 1000:1600;

clf; hold on;
plot(1:n*q, f1, 'LineWidth', 2);
plot(1:q:n*q,f, 'b.', 'MarkerSize', 25);
axis tight; box on;
axis([1000 1600 -0.09 0.045]);
set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
saveas(gcf, [rep 'zooming.eps'], 'epsc');


%%
% Same but after agressive sub-sampling.

r = 4; % subsampling factor
p = n/r;
g = f(1:r:end);
gh = fft(g);

gh1 = zeros(q*n,1); 
gh1(1:p/2) = gh(1:p/2);
gh1(end-p/2+2:end) = gh(p/2+2:end);
gh1(p/2+1) = gh(end/2+1)/2; 
gh1(end-p/2+1) = gh(end/2+1)/2;
g1 = q*r * real( ifft(gh1) ); % must be real anyway


clf; hold on;
plot(1:q:n*q,f, 'b.', 'MarkerSize', 20);
plot(1:q*r:n*q,g, 'r.', 'MarkerSize', 20);
plot(1:n*q, f1, 'b');
plot(1:n*q, g1, 'r');
axis tight; box on;


