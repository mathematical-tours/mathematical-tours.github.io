%%
% Test for TV denoinsing using a 4-fold finite differences gradient.
% Test for various regularization parameters lambda.

addpath('toolbox/');
addpath('imgs/');

% size of the image
N = 512/2;

if not(exist('name'))
name = 'lena';
name = 'square';
name = 'pacman';
name = 'square-tube-50';
name = 'cross';
end

rep = ['../results/total-variation/' name '/'];
if not(exist(rep))
    mkdir(rep);
end


%%
% Helpers and parameters.

dotp = @(x,y)sum(x(:).*y(:));
ampl = @(z)sqrt(sum(z.^2,3));
% saturation factor for display
rho = .05;

%%
% load the image

if strcmp(name, 'square-tube-50')
    P = round(N*.07);
    f0 = load_image(name, N-2*P);
    f = zeros(N); 
    f(P+1:N-P,P+1:N-P) = f0;
else
    f = load_image(name, N);
end

f = rescale(sum(f,3));
f = double(f>.5);
imwrite(f, [rep 'original.png'], 'png');

%%
% Test grad and div operator. div=-grad^T

% image and its Laplacian.
[grad,div,lapl] = load_grad(N);
figure(1);

clf;
imageplot({f lapl(f)});
% check for adjointness
e = abs( dotp(grad(f),grad(f)) + dotp(div(grad(f)),f) ) / abs(dotp(grad(f),grad(f)));
fprintf('Should be 0: %.1e.\n',e);
% check for Laplacian norm, should be 16
% [L,e] = compute_operator_norm(lapl,randn(N), 20);

%%
% Results with noise. Here no noise

sigma = 0;
randn('state', 1234);
y = f+randn(N)*sigma;
u = image_saturation(y,rho);
imwrite(u, [rep 'noisy-img.png']);


q = 10; % number of tests
% for N=256, [.1,.17]
lambda_list = linspace(.001,8,q); % 
options.niter = 500*20;

for i=1:length(lambda_list)
    lambda = lambda_list(i);
    [f1,E,z] = perform_tv_denoising(y,lambda, options);
    options.z0 = z;
    
    
    % display quantized colormap
    r = 14; % #levellines
    t = linspace(0,1,N);
    clf; hold on;
    imagesc(t,t,f1);
    contour(t,t,f1,linspace(min(f1(:)), max(f1(:)),r), 'k');
    colormap(parula(r-1));
    % caxis([0 1]);
    caxis([min(f1(:)) max(f1(:))]);
    axis image; axis off; axis ij; drawnow;
    saveas(gcf, [rep 'tv-flow-' num2str(i) '.png'], 'png');
end
