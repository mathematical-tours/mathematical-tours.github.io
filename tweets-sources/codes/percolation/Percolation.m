

addpath('../toolbox/');
rep = MkResRep();

n = 128*4;
% gaussian blur
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

M = GFilt(rand(n),2);
% track evolution of changes
tlist0 = linspace(.45,.8,500);
SZ = [];
for it=1:length(tlist0)
    progressbar(it,length(tlist0));
    I = M<tlist0(it);
    r = 5;
    I(end/2,end/2) = 1; % ensure not background
    L = bwlabel(I,4);
    SZ(it) = sum(L(:)==L(end/2,end/2));
end
plot(tlist0, SZ/n^2);


[SZ,J] = unique(SZ); tlist0 = tlist0(J);

q = 50;
tlist = interp1(SZ/n^2, tlist0, linspace(0.01,.99,q) );


for it=1:q
    s = (it-1)/(q-1);
    I = M<tlist(it);
    r = 5;
    % I = I(end/2-r-1:end/2+r,end/2-r-1:end/2+r);
    I(end/2,end/2) = 1; % ensure not background
    %
    L = bwlabel(I,4);
    J = find(L==L(end/2,end/2));
    %
    U = ones(n,n,3);    
    c = [(1-s), 0 s];
    for i=1:3
        V = ones(n);
        V(J) = c(i);
        U(:,:,i) = V;
    end
    %
    clf;
    imagesc(U); axis image; axis off; colormap gray;
    drawnow;
    imwrite(U, [rep 'anim-' znum2str(it,3) '.png'] );
end