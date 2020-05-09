% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it,name)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);

% initial points
m = 16*2;
if 1
    % on a circle
t = (0:m-1)'/m;
x0 = .9*exp(2i*pi*t);
else
x0 = (2*rand(m,1)-1) + 1i*(2*rand(m,1)-1);
x0 =linspace(-.8,.8,m)';
end

% time for integration
n = 2048*4;
t = linspace(-1,1,n)*100;

% sym matrix
Orth = @(th)[cos(th) sin(th); -sin(th) cos(th)];
SMat = @(th,eta)Orth(th)*diag([1 eta])*Orth(-th);



%  pure imaginary -> real positive
A0 = [0 1; -1 0]*.9;
A1 = -1.6*SMat(pi/3,.1);

% real positive -> real negative
A0 = 1.7*SMat(pi/3,.1); 
A1 = -1.7*SMat(pi/3+pi/4,.1);

% to display vector fields
ta = linspace(-1,1,12);
[Y,X] = meshgrid(ta,ta);

q = 70;
for it=1:q
    s = (it-1)/(q-1);
    A = (1-s)*A0 + s*A1;
    y = IntegrateLin(A,x0,t);
    % vector field
    V = reshape(A*[X(:)';Y(:)'], [2 length(ta) length(ta)]);
    V = permute(V, [2 3 1]);
    %
    clf; hold on;
    plot(real(y), imag(y), 'LineWidth', 1, 'Color', [s 0 1-s]);
    % plot(real(x0), imag(x0), 'k.', 'MarkerSize', 15);
    quiver(X,Y,V(:,:,1),V(:,:,2), 'k', 'filled', 'LineWidth', 1, 'AutoScaleFactor', 1.2);    
    axis equal;
    axis([-1 1 -1 1]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;  
    mysaveas(it,'anim');
end

% display evolution on the plane tr/det
D = []; T = [];
for it=1:q
    s = (it-1)/(q-1);
    A = (1-s)*A0 + s*A1;
    D(it) = det(A); T(it) = trace(A);
end
Tmax = 2;
T1 = linspace(-Tmax,Tmax,200);
D1 = linspace(-Tmax^2/4,Tmax^2/4,200);
for it=1:q
    s = (it-1)/(q-1);
    clf; hold on;
    plot(T,D, 'k', 'LineWidth', 2);
    plot(T1,T1.^2/4, 'k', 'LineWidth', 1);
    plot(T1,T1*0, 'k', 'LineWidth', 1);
    plot(D1*0,D1, 'k', 'LineWidth', 1);
    plot(T(it),D(it), '.', 'Color', [s 0 1-s], 'MarkerSize', 35);
    axis([-Tmax,Tmax,-Tmax^2/4,Tmax^2/4]); 
    axis off;
    drawnow;   
    mysaveas(it,'planar'); 
end


