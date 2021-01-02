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



% low pass
s = 1/2;
h{1} = zeros(n);
h{1}([end,1,2],[end,1,2]) = [0 s 0; s 1 s; 0 s 0] / (1+4*s);
% d/dx
h{2} = zeros(n);
h{2}([end,1,2],[end,1,2]) = [0 0 0; -1 0 1; 0 0 0]; 
% d/dy
h{3} = zeros(n);
h{3}([end,1,2],[end,1,2]) = [0 -1 0; 0 0 0; 0 1 0]; 
% laplacian
h{4} = zeros(n);
h{4}([end,1,2],[end,1,2]) = [0 -1 0; -1 4 -1; 0 -1 0]; 

resc = @(f).5+.5*f/max(abs(f(:)));

fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

imwrite( rescale(Upsc(f,8)), [rep 'input.png'] );

q = 70;
for it=1:q
    t = (it-1)/q;
    i = floor(t*3); s = (t*3)-i;
    h1 = h{i+1}*(1-s) + h{i+2}*s;
    %
    imwrite( resc(Upsc(h1([end,1,2],[end,1,2]),32)), [rep 'filter-' znum2str(it,2) '.png'] );
    F = {h1 h{1} h{1} h{1} h{1} h{1} h{1} h{1} h{1} h{1}}; 
    clf;
    g = f;
    for k=1:length(F)
        g = fconv(g,F{k});
        if k<7
        subplot(2,3,k);
        imagesc(g); 
        end
        if k==length(F)
            imwrite( rescale(Upsc(g,8)), [rep 'anim-' num2str(k) '-' znum2str(it,2) '.png'] );
        end
    end
    drawnow;
    
    
    %imwrite( resc(Upsc(fs,8)), [rep 'anim-' znum2str(it,2) '.png'] );
    %imwrite( resc(Upsc(fftshift(-g),8)), [rep 'filter-' znum2str(it,2) '.png'] );
end