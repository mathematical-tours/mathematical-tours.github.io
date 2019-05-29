%%
% Display DCT matrices

addpath('../toolbox/');
rep = MkResRep();

% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));


n = 32;
[Y,X] = meshgrid(0:n-1,0:n-1);

A1 = cos(X.*Y*pi/(n-1));
A1(:,[1,end]) = A1(:,[1,end])/sqrt(2);
A1([1,end],:) = A1([1,end],:)/sqrt(2);
imwrite(rescale(Upsc(A1)), [rep 'dct1.png']);

A2 = cos((X+1/2).*Y*pi/n);
A2(:,1) = A2(:,1)/sqrt(2);


A4 = cos((X+1/2).*(Y+1/2)*pi/n)*sqrt(2/n);

Cas= ( cos(X.*Y*2*pi/n) + sin(X.*Y*2*pi/n) ) * 1/sqrt(n);