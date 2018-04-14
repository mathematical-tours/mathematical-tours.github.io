%%
% 2D Metaballs

addpath('../toolbox/');
rep = MkResRep('2d');

n = 256;


K = 6; % #balls
% position of the centers

name = 'gauss';
name = 'inv';
name = 'invb';
name = 'inv2';
switch name
    case 'inv2'
        phi = @(d)1./(.01+d);
        T = 50;  % threshold for the balls, increase to decrease the size.
    case 'inv'
        phi = @(d)1./(.01+sqrt(d));
        T = 15;  % threshold for the balls, increase to decrease the size.
    case 'inva'
        phi = @(d)1./(.01+d.^1.1);
        T = -.15;  % threshold for the balls, increase to decrease the size.
    case 'invb'
        phi = @(d)1./(.01+d.^.3);
        T = -.15;  % threshold for the balls, increase to decrease the size.
    case 'gauss'
        s = .08;
        phi = @(d)exp(-d/(2*s^2));
        T = .3;  % threshold for the balls, increase to decrease the size.
end
T = -.15;

rand('seed', 1235);
w = rand(K,1)/2+1/2; % size
q = 100; % #movie steps

c0 = rand(K,1)+1i*rand(K,1); % center
s0 = 4*exp(2i*pi*rand(K,1)) .* (1+rand(K,1)) / 2; % linear speed
r0 = (rand(K,1)/2+1/2)*.5;

for i=1:q
    
    t = (i-1)/q;
    c = c0 + t*s0;   
    % c = c0 + exp(2i*pi*t*(-1).^(1:K)').*r0;
    
    
    RF = GenBall2D([real(c),imag(c)],w,n, phi, T);
    clf;
    imagesc(RF); 
    axis image; axis off;
    drawnow;
    imwrite(RF, [rep name '-' znum2str(i,3) '.png']);
end