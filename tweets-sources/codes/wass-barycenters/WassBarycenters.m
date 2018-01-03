addpath('./toolbox/');
addpath('./imgs/'); 

%% Read data

N = 200; % size of images
epsilon = 2; % for N=200

N = 100; % size of images
epsilon = .8; % for N=100


names = {'monge', 'kantorovitch'};
names = {'benamou', 'brenier'};

P = length(names);

vmin = .02;
vmin = .001;
H = zeros(N*N,P);
for i=1:P
    A = load_image(names{i}, N); A = sum(A,3);
    % saturate
    A = clamp(rescale(A),.2,.8);
    % inverse 
    A = rescale(-A,vmin,1);
    A = A / sum(A(:));
    H(:,i) = A(:);
end

rep = ['../results/wass-barycenters/' names{1}];
for i=2:P
    rep = [rep '-' names{i}];
end
rep = [rep '/'];
[~,~] = mkdir(rep);

n = size(H,1);
areaWeights = ones(n,1)/n;
H = H*n;

entropies = -sum(bsxfun(@times,H.*log(H),areaWeights),1);
maxEntropy = max(entropies);

%% Set up blur

imS = [N N];

if exist('imfilter') && 0
    % using image toolbox
    h = fspecial('gaussian',[1 max(imS)],epsilon);% hsize sigma
    h = h / sum(h);
    imBlur = @(x) imfilter(imfilter(x,h,'replicate'),h','replicate');
else
    blur = load_filtering('imgaussian', N);
    imBlur = @(x)blur(x,epsilon);
end

imagesc(imBlur(reshape(H(:,1),imS)));

blurColumn = @(x) reshape(imBlur(reshape(x,imS)),[],1);
blurAll2 = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));
blurAll = @(x) blurAll2(blurAll2(x));

% blurColumn = @(x) reshape(fastBlur(reshape(x,[targetSize,targetSize]),filterSize),[],1);
% blurAll = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));

imagesc(reshape(blurAll(H(:,1)),imS));
axis equal;
axis off;

%% weights


% bilinear interpolation
switch P
    case 4
        q=5;
        t = linspace(0,1,q);
        [T,S] = meshgrid(t,t); S = S(:); T = T(:);
        W = [(1-S).*(1-T) S.*(1-T) (1-S).*T S.*T]';
    case 2
        q = 20;
        t = linspace(0,1,q);
        W = [1-t; t];
end
Q = size(W,2);

%% Compute barycenter

entropyLimit = Inf; % no limit

steps = linspace(0,1,5);% linspace(-.5,1.5,9);

averages = cell(length(steps),length(steps));
barycenters = cell(length(steps),length(steps));

options.niter = 4000; %not 300
options.niter = 3000; %not 300
options.verb = 2;
options.tol = 1e-15;
resh = @(x)reshape(x, imS);
options.disp = @(x)imageplot(-resh(x));

col = [ [1 0 0]; [0 1 0]; [0 0 1]; [1 1 0]; [1/2 1 1/2]; [1/2 1/2 1] ];

B = {};
for i=1:Q
    w = W(:,i)'; w = w/sum(w);
    %
    B{i} = convolutionalBarycenter(H,w,areaWeights,blurAll,blurAll,entropyLimit,options);
    if isnan(sum(B{end}))
        warning('Problem');
        break;
    end
    U = rescale(-resh(B{i}));
    clf; imageplot(U); drawnow;
    % save image
    imwrite(U, [rep 'interp-' num2str(i) '.png'], 'png');
end

r = Q;
for i=Q-1:-1:2
    r = r+1;
    U = rescale(-resh(B{i}));
    clf; imageplot(U); drawnow;
    % save image
    imwrite(U, [rep 'interp-' num2str(r) '.png'], 'png');
end

