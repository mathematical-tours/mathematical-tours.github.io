%%
% Run k-means++

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


% generate mixtures
k0 = 30;
z0 = (.1+.1i) + .8*( rand(k0,1) + 1i*rand(k0,1) );

% mean/scale/anisotrop/orientation

p = 50; % #sample per cluster
gauss = @(m,s,a,t)m + s*exp(2i*t)*( randn(p,1) + 1i*randn(p,1)*a );

X = [];
a = .5; s = .03;
for k=1:k0
    X = [X; gauss(rand+rand*1i,s,a,rand)];
end
n = length(X);

clf;
plot(X, '.');
axis equal; axis([0,1,0,1]);
box on; set(gca, 'XTick', [], 'YTick', []);

mode = 'fp';
mode = 'keans++';

q = 30;
Col = distinguishable_colors(q);
Z = .5+.5i;
jt = 1; kt = 1;
for it=1:q
    D = abs(Z(:) - transpose(X(:)));
    [d,I] = min(D, [],1);
    % proba
    P = d.^2 / sum(d.^2);
    % farthest point
    switch mode
        case 'fp'
            % farthest point
            [~,i] = max(d);
        otherwise
            % kmeans++
            i = rand_discr(P, 1);
    end
    %
    if it==1
        Z = X(i);
    else
        Z(end+1) = X(i);
    end
    % coloring     
    D = abs(Z(:) - transpose(X(:)));
    [d,I] = min(D, [],1);
    P = d.^2 / sum(d.^2);
    %   
    s = ones(n,1)*12; % size
    clf; hold on;
    scatter( real(X), imag(X), s, Col(I,:), 'filled' );    
    for j=1:length(Z)
        plot(Z(j), '.', 'color', .5*Col(j,:), 'MarkerSize', 35);
    end
    axis equal; axis([0,1,0,1]);
    box on; set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    for m=1:2
        jt = jt+1;
        mysaveas('anim',jt)
    end
    % 
    m = P(:)/max(P);
    c = m*[1 0 0] + (1-m)*[0 0 1];
    clf; scatter( real(X), imag(X), s, c, 'filled' );
    axis equal; axis([0,1,0,1]);
    box on; set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    for m=1:2
        kt = kt+1;
        mysaveas('proba',kt)
    end    
end