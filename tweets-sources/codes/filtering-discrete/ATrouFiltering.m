%%
% Vizualization of a discrete filter.

% Translation-Invariant De-Noising
% R. R. CoifmanD. L. Donoho
% --> undecimated or stationary wavelet transform.

% alled algorithme à trous, introduced by Holschneider et al. [303]. 

% a trou algorithm 
% M. Holschneider, R. Kronland-Martinet, J.Morlet,and P. Tchamitchian.
% A real-time algorithm for signal analysis with the help of the wavelet transform. In Wavelets, Time-Frequency Methods and Phase Space, pp. 289?297. Springer-Verlag, 1989.

% M. J. Shensa. The discrete wavelet transform: Wedding the à trous and Mallat algorithms. IEEE Trans. Signal Process., 40(10):2464?2482, 1992.
% S. Mallat and S. Zhong. Characterization of signals from multiscale edges. IEEE Trans. Patt. Anal. Mach. Intell., 14(7):710?732, 1992.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)),:);

% binary image.
n = 256;
name = 'hibiscus';
f = load_image(name, n);
f = rescale(f);

% to display the PSF
if 0
    f = zeros(n);
    f(end/2,end/2) = 1; 
end

% low pass
s = 1/2;
% h = zeros(n);
% h([end,1,2],[end,1,2]) = 

h = [0 s 0; s 1 s; 0 s 0] / (1+4*s);
h = conv2(h,h);

% haar
h = ones(2);

% Wavelet daubechies 4
h = [.482962913145, .836516303738, .224143868042, -.129409522551]';

h = [.35355339059327, 0.70710678118655, 0.35355339059327 0]';

h = h*h';

% CDF 5/3

h = h/sum(h(:));

atrou = 1;

g = f;
H = 1;
jmax = 6;
for j=1:jmax
    % t = (it-1)/q;
    if atrou==1
        % h1 = zeros(size(h)*2^(j-1));
        % h1(2^(j-1):2^(j-1):end,2^(j-1):2^(j-1):end) = h;
        h1 = upsample(h,2^(j-1));
        
        imwrite( Upsc(rescale(-h1), n/size(h1,1) ), [rep 'filter-' znum2str(j,2) '.png'] );
        imwrite( Upsc(rescale(-H), n/size(h1,1) ), [rep 'globalfilt-' znum2str(j,2) '.png'] );
        
    else
        h1 = h;
    end
    clf;
    imagesc(g); axis image; axis off;
    colormap gray(256);
    drawnow;
    
    
    
    H = conv2(H,h1);
    
    
    g = imfilter(g,h1, 'replicate');
    if atrou==0
        g = g(1:2:end,1:2:end,:);
    end
        
    % imwrite( Upsc(g, n/size(g,1) ), [rep 'image-' znum2str(j,2) '.png'] );
    
end

cv = @(x,y)imfilter(x,y, 'replicate');
subs = @(x)x(1:2:end,1:2:end,:);
mynorm = @(x)norm(x(:));

f1 = subs(cv(subs(cv(f,h)),h));
f2 = subs(subs(cv(cv(f,h),upsample(h,2))));


cv = @(x,y)imfilter(x,y);
cv = @(x,y)imfilter(x,y, 'replicate');
mynorm( cv(subs(f),h) - subs(cv(f,upsample(h,2))) )
mynorm( subs( cv(subs(f),h) ) - subs( subs(cv(f,upsample(h,2))) ) )