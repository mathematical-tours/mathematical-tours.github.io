%%
% Iterative convolution / TCL.

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

n = 2048;
T = 2*[0:n/2, -n/2+1:-1]'/n;



name = 'sawtoothbox';
name = 'sawtooth';
name = 'boxes';
name = 'twoboxes';

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
a = mean(f0(T));
m = mean(f0(T).*T)/a;
s = sqrt(mean(f0(T).*(T+m).^2)/a);
s = 1;

f = @(tau)f0( (T/tau)*s+m );

q = 70; % frame anim
tlist = linspace(1,8,q);


for i=1:q
    T = tlist(i);
    s = (i-1)/(q-1);
    
    
    g = ifft( fft(f(1/sqrt(T))).^T );
    
    g = real(g);
    g = g/sum(g)*n;
    if 1 % i==1
        vmax = max(max(g));        
    end
    clf; hold on;
    plot(fftshift(g), 'LineWidth', 2, 'Color', [s 0 1-s]);
    axis tight; box on;
    axis([1 n -.15 vmax*1.03]);
    set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/3);
    saveas(gcf, [rep name '-' znum2str(i,2) '.png'], 'png');
    drawnow;
end
% 
