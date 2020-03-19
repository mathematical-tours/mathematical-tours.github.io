%%
% Just display gaussian convolution.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 2048;
q = 80;

% gaussian blur
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); % G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft( fft(x).*fft(y) ) );
GFilt = @(f,s)fconv(f, G(s));




f0 = zeros(n,1); p = 10;
f0(floor(rand(p,1)*n)) = sign(randn(p,1)) .* (6 + randn(p,1));
f0 = rescale(cumsum(f0), .05, .95 );
tlist = .01 + 900*linspace(0,1,q).^2;


f0 = .5 + 3*(rand(n,1)-1/2);
tlist = 5 + 300*linspace(0,1,q).^2;


x = linspace(0,1,n)';
mybox = @(m,s)abs(x-m)<s/2;
f0 = mybox(.3,.1) + .3*mybox(.5,.6) + 1.5*mybox(.8,.1);
f0 = rescale(f0, .05, .95 );
tlist = .01 + 600*linspace(0,1,q).^2;



f0 = sin(4*pi*x) + 4*mybox(.5,.5);
f0 = rescale(f0, .05, .95 );
tlist = .01 + 900*linspace(0,1,q).^2;


for it=1:q
    clf; hold on;
    for jt=it:it
        s = (jt-1)/(q-1);
        col = [s 0 1-s];
        f = GFilt(f0,tlist(jt));
        % h = plot(f, 'LineWidth', 2, 'Color', col);
        % alpha(h, 1/2);
        h = area(f, 'FaceColor', col, 'EdgeColor', 'none');
        h.FaceAlpha = 1;
    end    
    axis([1 n 0 1]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'FontSize', 20, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end