%%
% Approximation of signals using fft.

addpath('../toolbox/');
rep = MkResRep('signal');

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

n = 1024;
name = 'piece-regular';
f = rescale(load_signal(name, n));

m = [0:n/2,-n/2+1:-1]';
q = 70;
plist = round(linspace(1,n/4,q));

plist = round(1 + (n/2-1)*linspace(0,1,q).^2);

for i=1:q
    t = (i-1)/(q-1);
    p = plist(i);
    f1 = real(ifft( fft(f).*(abs(m)<=p) ) );
    clf;
    plot(1:n,f1,'LineWidth', 2, 'Color', [t 0 1-t]);
    axis([1 n -.1 1.02]);
    set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/2); drawnow;
    saveas(gcf, [rep name '-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, name '-');