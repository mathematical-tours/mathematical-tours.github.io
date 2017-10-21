%%
% Test for the fractional derivative of some input function.

rep = 'results/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
SetTickOff = @()set(gca, 'XTick',[], 'YTick',[]);


N = 2048*8*2;
t = [0:N/2,-N/2+1:-1]'/N;
x = linspace(-1,1,N)';

s = .01/4;
f = exp(-x.^2/(2*s^2));

kappa = @(a)(1i)^a * sign(t).*abs(t).^a;
kappa = @(a)(1i*t).^a;
fracD = @(f,a)real( ifft( fft(f) .* kappa(a) ) );

Q = @(x)x/max(abs(x));

a = 2;

alist = (0:15)*.2;

for i=1:length(alist)
    c=(i-1)/(length(alist)-1);
    a = alist(i);
    clf; hold on;
    plot(x,Q(fracD(f,a)), 'LineWidth', 2, 'Color', [c 0 1-c]);
    plot(x,x*0, 'k--');
    plot([0 0],[-1 1], 'k--');
    axis tight;
    axis([-8*s 8*s -1 1]);
    box on;
    SetAR(2/3);
    SetTickOff();
    legend(['\alpha=' num2str(a)]);
    set(gca, 'FontSize', 35);
    saveas(gcf, [rep 'fracder-' num2str(i) '.eps' ], 'epsc');
end


clf; hold on;
    plot(x,x*0, 'k--');
    plot([0 0],[-1 1], 'k--');
for i=1:length(alist)
    c=(i-1)/(length(alist)-1);
    a = alist(i);
    plot(x,Q(fracD(f,a)), 'LineWidth', 2, 'Color', [c 0 1-c]);
    axis tight;
    axis([-8*s 8*s -1 1]);
    box on;
    SetAR(2/3);
    SetTickOff();
    % legend(['\alpha=' num2str(a)]);
    set(gca, 'FontSize', 35);
end
    saveas(gcf, [rep 'fracder-all.eps' ], 'epsc');