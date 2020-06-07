
% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 1 + 2048; 
t = (0:n-1)'/n;

s = .06;

f = exp(-(t-1/2).^2/(2*s^2)) .* sin(12*pi*t);
f = exp(-(t-1/2).^2/(2*s^2)) .* (t-1/2);
f = exp(-(t-1/2).^2/(2*s^2)) .* (t-1/2).^2;
%
f = exp(-(t-1/2).^2/(2*s^2));


mypow = @(x,s)abs(x).^s .* exp(1i*sign(x)*s*pi/2);


if mod(n,2)==0
    om = [0:n/2, -n/2+1:-1]';
else
    om = [0:(n-1)/2, -(n-1)/2:-1]';
end
L = n^2*sin(pi/n*om).^2;
D = n*sin(pi/n*om);
Lapl = @(f,s)real( ifft( fft(f) .* L.^(s/2) ) );
Der = @(f,s)ifft( fft(f) .* mypow(D,s) );

q = 70; 
slist = linspace(1e-3,6,q);

%% 
% Laplacian.
operator = Der;
operator = Lapl;


F = []; 
mem = 30;
for it=1:q    
    r = (it-1)/(q-1);
    %
    fL = operator(f,slist(it)); 
    fL = fL-fL(1); fL = fL/max(abs(fL)); 
    F = [fL,F]; F = F(:,1:min(mem,end));
    % 
    clf; hold on;
    for j=size(F,2):-1:1
        m = .5+.5*(j-1)/(mem-1);
        plot(t, F(:,j), 'LineWidth', 1, 'color', [r 0 1-r]*(1-m) + [1 1 1]*m );
    end
    plot(t, F(:,j), 'LineWidth', 4, 'color', [r 0 1-r] );
	box on; set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'XTick', [], 'YTick', []);
    axis([.1 .9 -1.03 1.03]);
    drawnow;
    mysaveas(it);    
end


% replay

%%
% Deriv