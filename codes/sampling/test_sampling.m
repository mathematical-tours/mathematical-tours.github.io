%%
% test for Fourier transform, sampling, interpolation 

name = 'hibiscus';
F = imread([name '.png']);

F = sum(F,3)/255;
clf; imagesc(F); axis image; colormap gray;

% extract 1-D trace
f= F(:,end/2);
n = length(f);
clf; plot(f); axis tight;

% make it periodic 
t = linspace(0,1,n)';
f = f - f(1) + (f(1)-f(end))*t;
plot(f); axis tight;

% plot its fft
fh = fft(f);
plot( log(abs(fh)) ); axis tight;
plot( 2*pi*(-n/2:n/2-1)/n, log(abs(fftshift(fh))) ); axis tight;

% log/log plot for decay -> 1/n
clf;
plot(log(2:n/2),log(abs(fh(2:n/2))));

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

%%
% Compare with cubic splines


g_spline = interp1(1:q*r:n*q,g, 1:n*q, 'spline');

clf; hold on;
plot(1:q*r:n*q,g, 'r.', 'MarkerSize', 20);
plot(1:n*q, g1, 'r');
plot(1:n*q, g_spline, 'b');
axis tight; box on;

%%
% Display cubic spline kernel phi for dirac.
% WARNIG: the interpolation formula is *not* sum_n f(n) phi(x-n) !

meth = 'linear';
meth = 'pchip';
meth = 'nearest';
meth = 'spline';
m = 10;
a = (-m:m)'; b = zeros(2*m+1,1); b(m+1)=1;
ai =linspace(-m,m,1024*2); 
bi = interp1(a,b,ai, meth);
clf; hold on;
plot(ai,bi);
plot(a,b, 'r.', 'MarkerSize', 20);

% just b-spline

plot(x, eval_spline(3, x) );



