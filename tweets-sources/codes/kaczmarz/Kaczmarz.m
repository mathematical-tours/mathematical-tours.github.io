%%
% kaczmarz iterative projection vizualization.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name, it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


dotp = @(x,y)real(x.*conj(y));

% number of constraints 
k = 15;
% angle of the constraints
t = cumsum([rand(k,1)+.5]); t = -(pi*.5*rand/k) + (t/t(end)*pi);


n = 15; 
r = cumsum(2+rand(n,1)); r = rescale(r,.7,2.5); 
x = r .* exp( .03 * 2i*pi*randn(n,1) );

% x = 2*(randn(n,1) + 1i*randn(n,1));

q = 50;

col = distinguishable_colors(n);

for it=1:q    
    s = mod(it-1,k)+1; % #constraint
    % s = ceil(rand*k);
    clf; hold on;
    for i=1:k % display constraints
        plot([-2 2]*exp(1i*t(i)), 'k', 'LineWidth', 1);
    end
    plot([-2 2]*exp(1i*t(s)), 'k', 'LineWidth', 5);
    % projection on constraint
    x(:,end+1) = x(:,end) - dotp(x(:,end),1i*exp(1i*t(s))) * 1i*exp(1i*t(s));
    for j=1:n
        plot(x(j,:), '.-', 'LineWidth', 2, 'MarkerSize', 10, 'color', col(j,:));
        plot(x(j,end), '.', 'MarkerSize', 30, 'color', col(j,:));
    end
    axis equal; axis([-1 1 -1 1]); axis off;
    drawnow;
    % mysaveas('anim', it);
end