%%
% Texture synthesis by randomizatio of the phase.

name = 'wood';
name = 'wood1';
name = 'colors';

name = 'mnms';
name = 'fabric';
name = 'stone';

mynorm = @(x)norm(x(:));

addpath('../toolbox/');
rep = MkResRep(name);

n = 512;
f0 = rescale( load_image(name, n) );
f0 = Periodize(f0);
% f0 = clamp(f0);

w = randn(n)/n; 
w = w-mean(w(:))+1/n^2;
f = real(ifft2(fft2(f0).*repmat(fft2(w), [1 1 size(f0,3)])));

F = fft2(f0);
G = fft2(randn(n));
t = 0.5;

q = 40; 
for it=1:q
    t = (it-1)/(q-1);
    Theta = (1-t)*angle(F) + t*repmat(angle(G), [1 1 3]);
    f = abs(F) .* exp( 1i*Theta );
    f(1,1,:) = F(1,1,:);
    f = real(ifft2(f));
    imageplot(clamp(f));
    drawnow;
    imwrite(clamp(f), [rep 'anim-' znum2str(it,2) '.png'] );
end

