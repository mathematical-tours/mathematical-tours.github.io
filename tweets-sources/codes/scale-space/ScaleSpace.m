%%
% Progressive gaussian convolution.

n = 2048*2; 

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

t = [0:n/2,-n/2+1:-1]' /n;
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); 

% true periodization
t = (0:n-1)'/n;
G = @(s,m)exp(-(t-m).^2/(2*s^2)); 
G = @(s)normalize( G(s,0)+G(s,1) + G(s,-1)+G(s,2)+ G(s,-2)+G(s,3) );

fconv = @(x,y)real( ifft( fft(x).*fft(y) ) );
GFilt = @(f,s)fconv(f, G(s));


stdi = @(x)x/max(abs(x));
randn('state', 23);
x0 = randn(n,1); x0 = x0-mean(x0);

rho = 1.5;
q = 200;
smin = .01;
smax = .5;
slist = smin * (smax/smin) .^ linspace(0,1,q);

derorder = 1;
derorder = 2;
 
for it=1:q
    t = (it-1)/(q-1);
    s = slist(it);
    x = GFilt(x0,s); xi = stdi(x);
    clf; hold on;
    plot( xi, 'color', [t 0 1-t], 'LineWidth', 2  );
    % zero crossing of derivative
    switch derorder
        case 1
            D = x([2:end 1])-x;
            [u,v] = linzeros(1/2+(1:n), D, xi);
        case 2
            D = 2*x - x([2:end 1])-x([end 1:end-1]);
            [u,v] = linzeros((1:n), D, xi);
    end
    plot( u, v, '.', 'color', [t 0 1-t], 'MarkerSize', 30  );
    axis([1 n -1.03 1.03]); box on;  set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas('curves', it);
end

clf; hold on;
for it=1:q
    t = (it-1)/(q-1);
    s = slist(it);
    x = GFilt(x0,s); xi = stdi(x);
    % zero crossing of derivative
    switch derorder
        case 1
            D = x([2:end 1])-x;
            [u,v] = linzeros(1/2+(1:n), D, xi);
        case 2
            D = 2*x - x([2:end 1])-x([end 1:end-1]);
            [u,v] = linzeros((1:n), D, xi);
    end
    plot( u, it, '.', 'color', [t 0 1-t], 'MarkerSize', 15  );
    axis([1 n 0 q+1]); box on;  set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas('points', it);
end
