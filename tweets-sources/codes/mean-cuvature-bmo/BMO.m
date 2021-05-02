%%
% Bence?Merriman?Osher (BMO) algorithm
% Motion of multiple junctions: A level set approach
% Barry Merriman, James K Bence, Stanley Osher
% Diffusion generated motion by mean curvature; 1992
% B Merriman, J Bence, S Osher

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);

n = 512*2; 

% gaussian kernel
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt0 = @(f,s)fconv(f, G(s));
GFilt = @(f,s)imgaussfilt(f,s);
%

t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
f = double( max(abs(X),abs(Y))<.5 );

name = 'elephant';
name = 'triskel';
name = 'bunny';
f = load_image(name, n);
f = double( rescale(f,0,1)<.5 );

q = 100;
niter = 300;
disp_list = round( 1+(niter-1)*linspace(0,1,q).^3 );
disp_list = unique(disp_list);
idisp = 1;
s = 5*2;

for it=1:niter
    if it==disp_list(idisp)
    clf;
    imageplot(-f);
    drawnow;
    mysaveas(idisp);
    idisp = idisp+1;
    end
    f = GFilt(f,s);
    f = double( f>.5 );    
end