%%
% test for synthesis of dynamic Gaussian processes.

addpath('toolbox/');

rep = ['results/synthesis/'];
if not(exist(rep))
    mkdir(rep);
end

%%
% Helpers.

s = 2.5;
normalize = @(x)rescale( clamp( (x-mean(x(:)))/std(x(:)), -s,s) );
direction = @(rho,theta)rho*[cos(theta);sin(theta)];

%%
% Warm-up phase (in order to remove first frames).

WarmUp = 10;

%%
% Load kernel.

clear options;
n = 256; % space dimension.
options.nbframes = 120 + WarmUp; % time dimension.

%%
% Parameters.

SpacialType = 'isotropic-low';
SpacialType = 'natural-smooth';
SpacialType = 'isotropic';
SpacialType = 'natural-rough';
SpacialType = 'anisotropic';
%
MovementType = 'rotozoom';
MovementType = 'static';
MovementType = 'translation';
MovementType = 'dezoom';
MovementType = 'rotation';
MovementType = 'zoom';
% 
NonlinearityType = 'spiking';
NonlinearityType = 'triangle';
NonlinearityType = 'sawtooth';
NonlinearityType = 'binarize';
NonlinearityType = 'identity';

%%
% Select spacial kernel.

% default parameters : grating
options.theta = pi/3;
options.rho = 20;
options.srho = .1;
options.stheta = .1;
KerMode = 'gabor-angle';
switch SpacialType
    case 'anisotropic'
        options.stheta = .2;
        options.srho = 6;        
	case 'isotropic'
        options.stheta = 100;
        options.srho = 6;     
	case 'isotropic-low' % lower spacial freq
        options.stheta = 100;
        options.rho = 2;
        options.srho = 3;
    case 'natural-rough'
        KerMode = 'matern';
        options.alpha = 1;
        options.eta = .1;
    case 'natural-smooth'
        KerMode = 'matern';
        options.alpha = 3;
        options.eta = .1;
end
H = load_kernel(n, KerMode, options);


%%
% Time covariance

options.sigmat = 40;

%%
% Geometric movement

switch MovementType
	case 'static'
    case 'zoom'
        options.scale = 1/1.03;
    case 'rotozoom'
        options.scale = 1/1.03;
        options.rotation = .01;
    case 'dezoom'
        options.scale = 1.03;
    case 'translation'
        options.translation = [.1 -.05]*10;
    case 'rotation'
        options.rotation = .01;
end

%% 
% Perform synthesis. 

F = perform_ar_synthesis(H, options);

%%
% Non-Linearity

switch NonlinearityType
    case 'identity'
        M = @(x)x;
    case 'binarize'
        M = @(x)double(x>median(x(:)));
    case 'sawtooth'
        M = @(x)mod(x,.5);
    case 'triangle'
        M = @(x)max(1-abs( mod(x*3,2)-1 ), 0); 
    case 'spiking'
        M = @(x)max(1-abs( mod(x*3,2)-1 ), 0).^3; 
    otherwise
        error('Unknown non-linearity');
end    


%%
% Applies pointwise non-linearity.

F1 = normalize( M( normalize(F(:,:,WarmUp+1:end)) ) );

%%
% Save files.


rep = ['../results/motions-clouds/' SpacialType '-' MovementType '-' NonlinearityType '/'];
[~,~] = mkdir(rep);

for i=1:size(F1,3)
    imwrite(rescale(F1(:,:,i)), [rep MovementType '-' znum2str(i,3) '.png']);
end

