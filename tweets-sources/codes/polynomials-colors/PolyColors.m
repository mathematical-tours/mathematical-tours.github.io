%%
% display in colors of complex polynomials

addpath('../toolbox/');
rep = MkResRep('rational');

% grid
r = 1.5;
n = 501;
t = linspace(-r,r,n);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;

if not(exist('nbr'))
    nbr = 1;
end
save_inter = 1; % progressively save

P = []; Pc = []; Proots = []; Ppoles = [];

F = ones(n);
clf; hold on;
iter = 1;
while true
    if iter==1
        axis equal; axis([-r r -r r]);
    end
    [a,b,button] = ginput(1);
    switch button
        case 1
            Proots(end+1) = a+1i*b;
            F = F .* (Z-(a+1i*b));
        case 3
            Ppoles(end+1) = a+1i*b;
            F = F ./ (Z-(a+1i*b));            
        otherwise
            break
    end
    clf; hold on;
    PolyDisp(F,t);
    % roots plot/poles
    plot(real(Proots), imag(Proots), 'b.', 'MarkerSize', 30);
    plot(real(Ppoles), imag(Ppoles), 'r.', 'MarkerSize', 30);
    saveas(gcf, [rep num2str(nbr) '-' znum2str(iter,3) '.png']);    
    iter = iter+1;
end

nbr = nbr+1;


% AutoCrop(rep, [num2str(nbr) '-'])
% > convert -delay 1x3 1-*.png 1.gif % 3x slower

