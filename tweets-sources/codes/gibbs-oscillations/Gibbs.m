%%
% Display evolution of Gibbs oscilations.

n = 2048; 

name = 'box';
name = 'ramp';
name = 'sinejump';

addpath('../toolbox/');
rep = MkResRep(name);
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


% load signal
t = linspace(0,1,n);
switch name
    case 'sinejump'
        x = .3*sin(4*pi*t) + ( abs(t-.5)<=.2 );
    case 'ramp'
        x = (t-.5) .* ( abs(t-.5)<=.2 );
    case 'box'
        x = abs(t-.5)<=.2;
    case 'quadbump'
        x = (abs(t-.5)<=.3) .* ( 1 - 3*(t-.6).^2 );
end



x = rescale(x(:));
clf; 
plot(t,x);

q = 50;
m = 6;
rlist = linspace(.003.^(1/m),.2.^(1/m),q).^m;

for it=1:q
    r = rlist(it);
    s = (it-1)/(q-1);
    fr = [0:n/2,-n/2+1:-1]';    
    x1 = real(ifft( fft(x) .* ( abs(fr)<=(n/2*r) ) ));
    clf;
    plot(t, x1, 'LineWidth', 2, 'Color', [s 0 1-s]);
    axis([0 1 -.1 1.1]); drawnow;
    set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    mysaveas(it);
end
% AutoCrop(rep, 'anim-'); 