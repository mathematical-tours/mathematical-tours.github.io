%%
% Test for 1D OT using cumulative distributions.

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
SetTickOff = @()set(gca, 'XTick',[], 'YTick',[]);
SetTickOn = @()set(gca, 'XTick',[0 1/2 1], 'YTick',[0 1/2 1]);

N = 256;
t = linspace(0,1,N)';
gauss = @(m,s)exp(-(t-m).^2/(2*s).^2);
normalize = @(x)x/sum(x(:));

s = .6;
A = @(t)gauss(.3*s+t,.05*s) + .5*gauss(.6*s+t,.15*s);
B = @(t).5*gauss(1-s*.2-t,.04*s) + .8*gauss(1-s*.5-t,.05*s) + .5*gauss(1-s*.8-t,.04*s);
vmin = .025;
A = @(t)normalize(vmin + A(t));
B = @(t)normalize(vmin + B(t));

q = 50;

for it=1:q
    s = (it-1)/(q-1);
    a = A(s*.43);
    b = B(s*.43);
    
    plot(t,[a b]);
    
    clf;
    area(t, a, 'FaceColor', 'r', 'EdgeColor', 'r');
    axis tight; SetAR(1/2); SetTickOff();
    saveas(gca, [rep 'input-1-' znum2str(it,2) '.png']);
    
    clf;
    area(t, b, 'FaceColor', 'b', 'EdgeColor', 'b');
    axis tight; SetAR(1/2); SetTickOff();
    saveas(gca, [rep 'input-2-' znum2str(it,2) '.png']);
    
    % cumulative
    ca = cumsum(a);
    cb = cumsum(b);
    % inverse cumulatives
    ica = interp1(ca, t, t, 'spline');
    icb = interp1(cb, t, t, 'spline');
    % composition of function
    Tab = interp1(t, icb, ca, 'spline'); % icb o ca
    Tba = interp1(t, ica, cb, 'spline'); % ica o cb
    % should be close to Id
    Iaa = interp1(t, Tba, Tab, 'spline'); % Tba o Tab
    Ibb = interp1(t, Tab, Tba, 'spline'); % Tab o Tba
    
    lw = 2;
    clf; hold on;
    plot(t, ca, 'r', 'LineWidth', lw);
    plot(t, cb, 'b', 'LineWidth', lw);
    axis tight; box on; SetTickOff();  SetAR(1);
    % legend('C_{\mu}', 'C_{\nu}');
    set(gca, 'FontSize', 20);
    saveas(gca, [rep 'cumul-' znum2str(it,2) '.png'], 'png');
    
    clf; hold on;
    plot(t, ica, 'r', 'LineWidth', lw);
    plot(t, icb, 'b', 'LineWidth', lw);
    axis tight; box on; SetTickOff();  SetAR(1);
    % legend('C_{\mu}^{-1}', 'C_{\nu}^{-1}');
    set(gca, 'FontSize', 20);
    saveas(gca, [rep 'icumul-' znum2str(it,2) '.png']);
    
    clf; hold on;
    plot(t, Tab, 'color', [1/2 1/2 0], 'LineWidth', lw);
    plot(t, Tba, 'color', [0 1/2 1/2], 'LineWidth', lw);
    % legend('T', 'T^{-1}');
    axis tight; box on; SetTickOff();  SetAR(1);
    drawnow;
    % set(gca, 'FontSize', 20);
    saveas(gca, [rep 'transports-' znum2str(it,2) '.png']);
    
end

%% Barycentric interpolation 


% AutoCrop(rep, 'transports-')

