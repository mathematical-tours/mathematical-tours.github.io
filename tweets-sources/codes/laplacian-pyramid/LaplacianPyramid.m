%% 
% Laplacian pyramid.

addpath('../toolbox/');
rep = MkResRep();

S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));


name = 'hibiscus';

n = 512/4; 
f = rescale(sum(load_image(name,n),3));


% gaussian blur
t = @(n)[0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s,n)normalize(exp(-t(n).^2/(2*s^2))); G = @(s,n)G(s,n)*G(s,n)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s,size(f,1)));
s = .8; % scale for the filtering

imwrite(f, [rep 'original.png']);
    
% Generate the pyramid
J = 5;
fj = f; 
D = {}; F = {};
for j=1:5
    % filtering
    fjOld = fj;
    fj0 = GFilt(fj,s);
    % sub-sampling
    fj = fj0(1:2:end,1:2:end);
    fjU = upsample(fj);
    %
    D{j} = fjOld-fjU;
    F{j} = fj;
	% display
    imageplot({F{j} D{j}});
    % save images
    imwrite(Upsc(fj,2^j), [rep 'gausspyr-' num2str(j) '.png']);
    m = 4;
    g = clamp(D{j}/std(D{j}(:)),-m,m);
    imwrite(Upsc(rescale(g),2^(j-1)), [rep 'laplpyr-' num2str(j) '.png']);
end
D{end+1} = fj;

% invert
fj = D{end};
for j=5:-1:1
    fjU = upsample(fj);
    fj = D{j}+fjU;    
end