%%
% Display Walsh matrix. 

addpath('../toolbox/');
rep = MkResRep();


W0 = [1 -1;1 1]/sqrt(2);

% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));


n = 8;
W = W0;
for i=1:p
    Wa = Upsc(W,2^n / 2^i);
    imwrite(rescale(Wa), [rep 'walsh-' znum2str(2^i,3) '.png']);
    clf;
    imagesc(Wa); axis image; axis off;
    drawnow;
    W = kron(W,W0);
end