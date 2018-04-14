%%
% Vornoi for Lp cost.

addpath('../toolbox/');
rep = MkResRep();

% second "continuous measure"
n = 500; % size of the image
N = n*n;

name = 'gaussian';
name = 'twobumps';

t = linspace(0,1,n);
[V,U] = meshgrid(t,t);
switch name
    case 'gaussian'
        z = [.7 .6];
        s = .12;
        mu = .02 + exp( (-(U-z(1)).^2-(V-z(2)).^2)/(2*s^2) );
        mu = mu/sum(mu(:));
    case 'twobumps'
        Z = {[.6 .9] [.4 .1]};
        S = [.05 .07];
        mu = ones(n)*.01;
        for k=1:length(Z)
            z = Z{k}; 
            s = S(k);
            mu = mu  + exp( (-(U-z(1)).^2-(V-z(2)).^2)/(2*s^2) );
        end
        mu = mu/sum(mu(:));
end

% first measure mu
if not(exist('Y0'))
    % click and play
    Y = [];
    clf; hold on;
    imagesc(t,t,-mu'); axis equal;
    while true
        axis([0 1 0 1]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 20);
        if button==3
            break;
        end
        Y(:,end+1) = [a;b];
    end
end
k = size(Y,2);
nu = ones(k,1)/k;

% coloring
C = distinguishable_colors(k+1);
C(4,:) = [];





clf; hold on;
imagesc(t,t,-mu);
s = 60*ones(k,1); % size
scatter( Y(2,:), Y(1,:), s, .8*C, 'filled' );
axis equal; axis([0 1 0 1]); axis off;
colormap gray(256);
saveas(gcf, [rep name '-measures.png'], 'png');


% gradient descent
tau = .05; % step size
%
niter = 200;
ndisp = 2;
plist = linspace(.3,5,niter); 
w = zeros(k,1);
for i=1:niter
    J = PowerDiagram(Y,n,w,1+(w-min(w))*15);
    % sum area captured
    for j=1:k
        r(j) = sum(sum( (J==j).*mu ));
    end
    % gradient ascent
    w = w - tau*( nu-r(:) );
    drawnow;
    if mod(i,ndisp)==1
        saveas(gcf, [rep name '-' znum2str((i-1)/ndisp+1,3) '.png'], 'png');
    end
end

% AutoCrop(rep, [name '-']); 

% convert iter-*.png lloyd.gif