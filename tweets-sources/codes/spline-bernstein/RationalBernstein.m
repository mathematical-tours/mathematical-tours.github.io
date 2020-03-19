
addpath('../toolbox/');
rep = MkResRep();


clf; hold on;
Z = [];
while true
    axis equal; axis([0 1 0 1]);  % axis off;
    [a,b,button] = ginput(1);  
    if button==3
        break;
    end
    plot(a,b, '.', 'MarkerSize', 25); 
    Z(end+1) = a + 1i*b;
    plot(Z, 'r', 'LineWidth', 2);
end
n = size(Z,2);
Z = Z(:);


% w = abs(linspace(-1,1,n)).^;
wmax = 40;
wlist = 1+linspace(0,1 ,q).^3*(wmax-1);

q = 50;
for it=1:q
    s = (it-1)/(q-1);
    %
    w = ones(n,1);
    w([3 end-3]) = [1 1]*wlist(it);
    % the basis
    B = ( BernsteinBasis(x,n) .* w' ) ./ ( BernsteinBasis(x,n) * w );
    clf;
    plot(x, B, 'LineWidth', 2); 
    axis([0 1 0 1]);
    box on;set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
    saveas(gcf, [rep 'basis-' znum2str(it,2) '.png']);
    drawnow;
    %
    y = ( BernsteinBasis(x,n) * (w.*Z) ) ./ ( BernsteinBasis(x,n) * w );
    %
    clf; hold on;
    plot(Z, 'k', 'LineWidth', 1, 'MarkerSize', 20);
    sz = 30 + w/max(w)*50; % size
    c = (w-1)/(wmax-1); c =  c * [1 0 0] + (1-c)*[0 0 1];
    scatter( real(Z), imag(Z), sz, c, 'filled' );
    plot(y, 'r', 'Color', [s 0 1-s], 'LineWidth', 2);
    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    drawnow;
end


