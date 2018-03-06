%%
% Test for 2D reaction diffusion.
% https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system
% https://en.wikipedia.org/wiki/The_Chemical_Basis_of_Morphogenesis
% https://blogs.mathworks.com/graphics/2015/03/16/how-the-tiger-got-its-stripes/
%  Gray-Scott model. http://mrob.com/pub/comp/xmorphia/
% http://www.karlsims.com/rd.html

if not(exist('name'))
    name = 'elongated';
end

rep = ['../results/reaction-diffusion/new/'  name '/'];
[~,~] = mkdir(rep);
addpath('../toolbox/');


Delta = @(f)-f ...
    + .20*(circshift(f,[ 1, 0]) + circshift(f,[-1, 0])  ...
    +      circshift(f,[ 0, 1]) + circshift(f,[ 0,-1])) ...
    + .05*(circshift(f,[ 1, 1]) + circshift(f,[-1, 1])  ...
    +      circshift(f,[-1,-1]) + circshift(f,[ 1,-1]));

%% GOOD PARAM %%

% Diffusion rates
da = 1;
db = .5;

switch name
    case 'elongated'
        f=.055;
        k=.062;
        tmax = 8000;
    case 'blobs'
        f=.055*.5;
        k=.062;
        tmax = 8000;
    case 'mixed'
        f=.055*.6;
        k=.062;
        tmax = 8000;
    case 'panther'
        f = 0.026;
        k = 0.053;
        tmax = 3000;
    case 'thin'
        f = 0.055;
        k = 0.063;
        tmax = 15000;
end


% Size of grid
n = 128;
n = 128*2;

% 5,000 simulation seconds with 4 steps per simulated second
dt = .25;
niter = ceil(tmax/dt);
ndisp = 80; 
kdisp = ceil(niter/ndisp);

% Initialize A to one
A = ones(n);
m = 8; % #seeds
rand('state', 12356);
I = randperm(n*n); 
B = zeros(n);
B(I(1:m))=1;

B = zeros(n);
B(end/2-1:end/2+1, end/2-1:end/2+1) = rand(3);
q = 1; 
B(end/2-q:end/2+q, end/2-q:end/2+q) = rand(2*q+1)>.5;


t = 0;
nframes = 1;
for i=1:niter
    [A,B] = deal( ...
        A + (da*Delta(A) - A.*B.^2 + f*(1-A))*dt, ...
        B + (db*Delta(B) + A.*B.^2 - (k+f)*B)*dt );
    % t = t+dt;
    if mod(i,kdisp)==1
        clf; imagesc(B);
        axis image; axis off; colormap parula(256)
        drawnow;
        % save image do disk
        B1 = apply_colormap(B,parula(256));
        imwrite(B1, [rep num2str( (i-1)/kdisp + 1 ) '.png']);
    end
end