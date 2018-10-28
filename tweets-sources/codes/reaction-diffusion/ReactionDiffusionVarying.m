%%
% Test for 2D reaction diffusion, using non-homogeneous weights.


name = 'bunny';
name = 'cat';
name = 'elephant';
name = 'france';

addpath('../toolbox/');
rep = MkResRep(name);


Delta = @(f)-f ...
    + .20*(circshift(f,[ 1, 0]) + circshift(f,[-1, 0])  ...
    +      circshift(f,[ 0, 1]) + circshift(f,[ 0,-1])) ...
    + .05*(circshift(f,[ 1, 1]) + circshift(f,[-1, 1])  ...
    +      circshift(f,[-1,-1]) + circshift(f,[ 1,-1]));

%% GOOD PARAM %%

% Diffusion rates
da = 1;
db = .5;

% elongated / blobs

%f= [0.026,0.026];
%k = [0.053,0.053];
tmax = 8000;
tmax = 7000;

switch name
    case 'cat'
        f = [0.026, .055*.7];
        k = [0.053, .062];
    case 'elephant'
        f = [0.026, .055*.5];
        k = [ 0.053, .062];
        %
        f = [.055, .055*.5];
        k = [.063, .062];
        tmax = 12000;
    case 'bunny'
        f = [.055, 0.026];
        k = [.063, 0.053];
        %
        f = [.055, 0.029];
        k = [.063, 0.055];
        tmax = 10000;
    case 'france'
        f = [.055, .055*.5];
        k = [.063, .062];
        tmax = 16000;


% Diffusion rates
da = 1/3;
db = .5/3;
end




% Size of grid
n = 128;
n = 128*2;

W = zeros(n,n);
W(1:end/2,:) = 1;



W = double(rescale(sum(load_image(name,n),3))>.5);

F = (1-W)*f(1) + W*f(2);
K = (1-W)*k(1) + W*k(2);

% 5,000 simulation seconds with 4 steps per simulated second
dt = .25;
niter = ceil(tmax/dt);
ndisp = 60;
kdisp = ceil(niter/ndisp);

% Initialize A to one
A = ones(n);
m = 4; % #seeds
rand('state', 12356);
I = randperm(n*n);
B = zeros(n);
B(I(1:m))=1;


t = 0;
nframes = 1;
for i=1:niter
    [A,B] = deal( ...
        A + (da*Delta(A) - A.*B.^2 + F .* (1-A))*dt, ...
        B + (db*Delta(B) + A.*B.^2 - (K+F).*B)*dt );
    % t = t+dt;
    if mod(i,kdisp)==1
        clf; imagesc(B);
        axis image; axis off; colormap parula(256)
        drawnow;
        % save image do disk
        if 1
        B1 = apply_colormap(B,parula(256));
        imwrite(B1, [rep znum2str( (i-1)/kdisp + 1,2 ) '.png']);
        end
    end
end
