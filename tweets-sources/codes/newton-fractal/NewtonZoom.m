%%
% Displays the bassin of attraction of Newton method on the complex plane.

addpath('../toolbox/');
rep = MkResRep('zoom');

if not(exist('q'))
    q = 5; % order of the polynomial
end
if not(exist('theta'))
    % basic newton
    rho = 0;
    %
    theta = 1.4; rho = .7;
end

f  = @(z)z.^q-1;
Df = @(z)q*z.^(q-1);
a = 1.5;
R = exp(2i*pi/q*(0:q-1)');

tau = 0; % avergaging
Newton = @(z)z - ( 1 - rho*exp(1i*theta) ) * f(z)./Df(z);
Phi = @(z)(1-tau)*Newton(z)+tau*z;

% grid size
N = 512;

% zoom point
p0 = exp(2i*pi/(2*q));
u0 = .393 + 1i*.762; % target center, q=3
u0 = .358 + 1i*.577; % target center, q=5
zf = 1.1; % zoom factor per image



% ininitial width of the window
A = a;
% initial center
u = 0; 

nimg = 50;  % #images
for iz = 1:nimg
    tz = (iz-1)/(nimg-1);    
    % 
    % update window location
    A = A/zf;
    u = u/zf + u0*(1-1/zf);
    % generate image
    [C1,x,y] = PlotNewtonFractal(Phi,R,N,A,u);
    % imwrite(C1, [rep 'newton-fractal-' num2str(q) '.png'], 'png');
    imwrite(C1, [rep 'newton-' num2str(q) '-' znum2str(iz,2) '.png'], 'png');
    %
    clf; hold on;
    imagesc(x,y,permute(C1, [2 1 3]));
    axis image; axis off;
    drawnow;
end
% saveas(gcf,  [rep 'newton-fractal-' num2str(q) '-pt.png'], 'png');