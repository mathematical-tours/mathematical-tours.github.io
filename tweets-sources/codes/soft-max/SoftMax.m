%%
% Plot of soft-max

name = 'sparsemax';
name = 'softmax';

addpath('../toolbox/');
rep = MkResRep(name);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);
rescale = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));

% stabilized log-sum-exp and soft max
max2 = @(S)repmat(max(S,[],2), [1 size(S,2)]);
LSE = @(S)log( sum(exp(S), 2) );
LSE = @(S)LSE( S-max2(S) ) + max(S,[],2);
SM = @(S)exp(S) ./ repmat( sum(exp(S),2), [1 size(S,2)]);
SM = @(S)SM(S-max2(S));

% sparse max
tsparse = @(x,y)clamp((x-y+1)/2, 0,1);
SparseMax = @(x,y)deal( tsparse(x,y), 1-tsparse(x,y)  );

% Project on sum(A)=A'*1=1.
ProjSimplex = @(y)max(bsxfun(@minus,y,max(bsxfun(@rdivide,cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0);

q = 200;
t = linspace(0,1,q);
[Y,X] = meshgrid(t,t);
X = X(:); Y = Y(:);
U = [X(:), Y(:)];

%% In 2D, display the max
eps_list = [1 .5 .1 .01 .001];
Q = 50;
eps_list = linspace(.001, .5,Q);
vmax = 1.2;
for i=1:Q
    c = (i-1)/(Q-1);
    epsilon = eps_list(i);
    switch name
        case 'softmax'
            R = LSE(U/epsilon)*epsilon;
        case 'sparsemax'
            % Sparse Max
            [X1,Y1] = SparseMax(X/epsilon,Y/epsilon);
            % Value of the max
            R = X.*X1 + Y.*Y1 - epsilon/2*( X1.^2 + Y1.^2 ) + epsilon/2;      
    end
    R = reshape(R,[q q]);
    % display
    r = 15; % #levellines
    clf; hold on;
    imagesc(t,t,R');
    contour(t,t,R',linspace(vmin,vmax,r), 'k');
    colormap(parula(r-1));
    caxis([vmin vmax]);
    plot([0 0 1 1 0 0], [0 1 1 0 0 1], 'color', [c 0 1-c], 'LineWidth', 5);
    axis image; axis off;
    saveas(gcf, [rep '2d-' znum2str(i,2) '.png']);
end

%% for a single histogram, display the argmax
n = 12;
randn('state', 1234);
a = rescale(rand(1,n))*.5+.2; a([1 3 6 7 11]) = [.97 1 .9 .93 .95];
for i=1:Q
    epsilon = eps_list(i);    
    switch name
        case 'softmax'
            b = SM(a/epsilon);
        case 'sparsemax'
            b =   ProjSimplex(a'/epsilon)';
    end
    % display
    c = (i-1)/(Q-1);
    clf; hold on;
    bar(1:n,b, 'FaceColor', [c 0 1-c]);
    plot(rescale(a), 'k.--', 'MarkerSize', 25, 'LineWidth', 1);
    axis([.5 n+.5 0 1]);
    box on; set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/3); box on;
    % saveas(gcf, [rep 'sm-' num2str(epsilon) '.eps'], 'epsc');
    saveas(gcf, [rep 'histo-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, ['2d-']); 
% AutoCrop(rep, ['histo-']); 

