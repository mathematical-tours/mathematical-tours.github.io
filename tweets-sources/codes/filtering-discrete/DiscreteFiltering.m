%%
% Vizualization of a discrete filter.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));


% binary image.

n = 32;
name = 'bunny';
f = load_image(name, n);
f = rescale(sum(f,3));
f = f>.5;


% gaussian blur
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));

G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
smin = 1e-3;
smax = 4;
resc = @(u)rescale(u);

if 1
Lapl = @(h)4*h - h([2:end,1],:) - h([end,1:end-1],:) - h(:,[2:end,1]) - h(:,[end,1:end-1]); 
G = @(s)Lapl(G(s));
smin = .01;
smax = 3;
resc = @(f).5+.5*f/max(abs(f(:)));
end


fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

imwrite( rescale(Upsc(f,8)), [rep 'input.png'] );

q = 50;
slist = smin + (smax-smin)*linspace(0,1,q);
for it=1:q
    s = slist(it);
    g = G(s);
    fs = fconv(f, g);
    fs = fs/max(abs(fs(:)));
    clf; 
    imagesc(fs);
    % imagesc(fftshift(g))
    % caxis([-1 1]);
    drawnow;    
    
    imwrite( resc(Upsc(fs,8)), [rep 'anim-' znum2str(it,2) '.png'] );
    imwrite( resc(Upsc(fftshift(-g),8)), [rep 'filter-' znum2str(it,2) '.png'] );
end


