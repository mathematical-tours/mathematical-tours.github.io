% Do robust least squares.


addpath('../toolbox/');
rep = MkResRep();


% #points
n = 15*4;
% #outliers
m = 20;

% #points
n = 15*2;
% #outliers
m = 8;

% data, well spread
t = cumsum( .5 + rand(n,1) );
t = rescale(t, .05, .95);


ymax = 1;
y = 1.5*(t-1/2);


%
ymax = .05;
y = (t-1/2) .* (t-.9) .* (t-.1) + ymax * .0 * randn(n,1) ;

% add a few outlier
I = 2 + randperm(n-4); I = I(1:m);
I = round(linspace(3,n-4,m));
y(I) = (-1).^(1:m)' .* (.7 + .3*rand(m,1))*ymax;


% degree for fitting
d = 14;



% features
x = t .^( 0:d );
%
% for display
T = linspace(0,1,512)';
X = T .^( 0:d );


% RW least squares



niter = 1000;
q = 50; % display
plist = linspace(2.5,1,q);

eta = ones(n,1);
w = zeros(d+1,1);
epsilon = 1e-5;
tau = 1;

W = [];
for it=1:q
    p = plist(it);
    % min |
    for i=1:niter
        eta1 = (epsilon +  abs(y-x*w)).^(p-2);
        eta = (1-tau)*eta + tau*eta1;
        % min < diag(eta)*(y-X*w), y-X*w>
        % (X'*diag(eta)*X) * w = X' * diag(eta)*y
        w = pinv( x'*diag(eta)*x ) * ( x' * diag(eta)*y );
    end  
    W(:,end+1) = w;
    %
    clf; hold on;
    for j=1:2:it
        s = (j-1)/(q-1);
        plot(T, X*W(:,j), 'color', .3*[s 0 1-s]+.7, 'LineWidth', 1);
    end
    plot(T, X*W(:,it), 'color', [s 0 1-s], 'LineWidth', 3);
    plot(t, y, '.k', 'MarkerSize', 15);
    box on;
    set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'XTick', [], 'YTick', []);
    axis([0 1 -ymax ymax]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end