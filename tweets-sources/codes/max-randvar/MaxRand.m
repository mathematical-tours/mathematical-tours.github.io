%%
% Display the evolution of the PDF of the max of n random variables.

addpath('../toolbox/');
rep = MkResRep();

n = 1025;  % for display

name = 'unif';
name = 'gauss';

a = 3; 
t = linspace(-a,a,n);

% cdf
switch name
    case 'unif'
        b = 2; 
        F = @(x) (x+b)/(2*b) .* (abs(x)<=b) + (x>b);
    case 'gauss'
        b = 2; 
        F = @(x) erf(x)/2+1/2;
end

% derivative
deriv = @(u)(u - u([1 1:end-1]))*n/(2*a);

q = 70; 
nmax = 20;
for i=1:q
    s = (i-1)/(q-1);
    F1 = deriv( F(t).^( 1 + (nmax-1)*s^2 ) );
    F1 = F1 / ( (2*a/n) * sum(F1) );
    clf; c = [s 0 1-s];
    area(t, F1, 'EdgeColor', c, 'FaceColor', c);
    axis([-a a 0 2]);
    set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 2/3 1]);
    drawnow;
    saveas(gcf, [rep name '-' znum2str(i,2) '.png']);
end