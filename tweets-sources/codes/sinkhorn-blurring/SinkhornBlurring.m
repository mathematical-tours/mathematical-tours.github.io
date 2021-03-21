
addpath('../toolbox/');
addpath('./mexEMD/');
rep = MkResRep();


nX = 12;
nY = nX;

X = exp(2i*pi/3) * ( randn(nX,1) + 1i*randn(nX,1)*.3 ) - 2;
Y = exp(2i*pi/3) * ( randn(nY,1) + 1i*randn(nY,1)*.3 ) + 2;
myax = [-3, 3, -1.8, 1.8]*1.7;
% myax = 2*[min(real([X;Y])), max(real([X;Y])), min(imag([X;Y])), max(imag([X;Y]))];

clf; hold on;
plot(X, 'r.', 'MarkerSize', 20);
plot(Y, 'b.', 'MarkerSize', 20);
axis equal; 
box on; set(gca, 'XTick', [], 'YTick', []);
axis(myax);

q = 50;

elist = ( .01 + .5*linspace(0,1,q).^2 ).^2; elist(1)=0;
for it=1:q
    epsilon = elist(it);
    [I,J,w] = myot(X,Y,epsilon);
    clf; hold on;
    for k=1:length(I)
        i = I(k); j = J(k);
        plot([X(i),Y(j)], 'color', .5*w(k)*[1 1 1], 'LineWidth', 1);
    end
    plot(X, 'r.', 'MarkerSize', 20);
    plot(Y, 'b.', 'MarkerSize', 20);
    axis equal;
    box on; set(gca, 'XTick', [], 'YTick', []);
    axis(myax);
    drawnow;
    saveas(gcf, [rep  'sinkhorn-' znum2str(it,2) '.png'], 'png');
end

% now with blurring
m = 20;
WX = randn(nX,m) + 1i*randn(nX,m);
WY = randn(nY,m) + 1i*randn(nY,m);
smax = 1;
for it=1:q
    s = smax*(it-1)/(q-1); % noise deviation
    X1 = X + s*WX; X1 = X1(:);
    Y1 = Y + s*WY; Y1 = Y1(:);
    [I,J,w] = myot(X1,Y1,0);
    clf; hold on;
    for k=1:length(I)
        i = I(k); j = J(k);
        plot([X1(i),Y1(j)], 'color', .5*[1 1 1], 'LineWidth', 1);
    end
    %plot(X, 'r.', 'MarkerSize', 20);
    %plot(Y, 'b.', 'MarkerSize', 20);
    plot(X1, 'r.', 'MarkerSize', 10);
    plot(Y1, 'b.', 'MarkerSize', 10);
    axis equal;
    box on; set(gca, 'XTick', [], 'YTick', []);
    axis(myax);
    drawnow;
    saveas(gcf, [rep  'blurring-' znum2str(it,2) '.png'], 'png');
end
    