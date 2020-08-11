%%
% Display the evolution of the 2D median

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


if not(exist('k'))
    k = 4; % #point
end
if not(exist('m'))
    m = 30; % #iteration
end

randn('state', 1223); rand('state', 123);
x = rand(k,1) + 1i*rand(k,1);
[~,I] = sort(angle(x-.5-.5i));
x = x(I);

% regular gon
x = .5+.5i + .4*exp( (0:k-1)'/k*2i*pi );


q = 150;
eta = .02; 
v = randn(k,1) + 1i*randn(k,1); 
v = eta*v./abs(v);
for it=1:q    
    clf; hold on;    
    y = x; 
    for j=1:m
        s = j/m;
        y = (y+y([2:end 1]))/2;
        a = 1;
        plot(y([1:end 1]), '-', 'color', [s^a 0 1-s^a], 'LineWidth', 1);
    end
    plot(x([1:end 1]), 'b.-', 'MarkerSize', 20, 'LineWidth', 2);
    plot(mean(x), 'r.', 'MarkerSize', 25);
    axis equal; axis([-.03 1.03 -.03 1.03]); 
    box on;  set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    % animate outliers
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
    mysaveas(['a-' num2str(k) '-'], it);
end