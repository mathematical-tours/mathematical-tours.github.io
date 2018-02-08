%%
% Iterative convolution / TCL.


rep = '../results/tcl/';
[~,~] = mkdir(rep);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

n = 2048;

t = 2*[0:n/2, -n/2+1:-1]'/n;


name = 'twoboxes';
name = 'boxes';
name = 'sawtooth';
name = 'sawtoothbox';


% box
bb = @(x,c,w)double(abs(x-c)<=w);
% sawtooth
st = @(x,c,w)(abs(x-c)<w) .* (x-c+w);


switch name
    case 'twoboxes'
        f0 = @(x)1.5*bb(x,-.3,.15) + bb(x,.3,.2);
    case 'boxes'
        f0 = @(x)bb(x,0,1/2);
    case 'sawtooth'
        r = -.1;
        f0 = @(x)st(x,.2,.2) + st(-x,.2,.2);
    case 'sawtoothbox'
        f0 = @(x)st(x,.3,.15)+.2*bb(x,-.3,.2);
end

% centers
a = mean(f0(t));
m = mean(f0(t).*t)/a;
s = sqrt(mean(f0(t).*(t+m).^2)/a);
s = 1;

f = @(tau)f0( (t/tau)*s+m );

K = 7;
clf; hold on;
vmax = -Inf;
for k=1:K
    s = (k-1)/K;
    g = ifft2( fft(f(1/sqrt(k))).^k );
    g = real(g);
    g = g/sum(g)*n;
    vmax = max(vmax,max(g));
    plot(fftshift(g), 'LineWidth', 2, 'Color', [s 0 1-s]);
end
axis tight; box on;
axis([1 n -.05 vmax*1.03]);
set(gca, 'XTick', [], 'YTick', []);
SetAR(1/3);
saveas(gcf, [rep name '.eps'], 'epsc');


clf; hold on;
vmax = -Inf;
for k=1:K
    s = (k-1)/K;
    g = ifft2( fft(f(1/sqrt(k))).^k );
    g = real(g);
    g = g/sum(g)*n;
    vmax = max(vmax,max(g));
    u = fftshift(g);
    U = cumsum(u); U = U/U(end);
    plot(U, 'LineWidth', 2, 'Color', [s 0 1-s]);
end
axis tight; box on;
axis([1 n -.03 1.03]);
set(gca, 'XTick', [], 'YTick', []);
SetAR(1/3);
saveas(gcf, [rep name '-repart.eps'], 'epsc');