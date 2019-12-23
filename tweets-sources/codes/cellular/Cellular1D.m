%%
% 1D celullar automata

addpath('../toolbox/');
rep = MkResRep();


% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));

% convert to binary 
tobin = @(x)double(dec2bin(x))-48;

% rule #
r = 126;
r = 30; 
R = tobin(r)';
R = [zeros(8-length(R),1);R];
R = R(end:-1:1);

n = 128/2; % size
k = n*2; % number of iterations


% initialization
X = zeros(n,k);
% X(round([.3 .5 .8]*n),1) = 1;
X(n/2,1) = 1;
X(:,1) = randn(n,1)>0;

imwrite(Upsc(rescale(1-X'),2), [rep 'anim-001.png']);
for it=2:k
    c = X([end 1:end-1],it-1) + 2*X(:,it-1) + 4*X([2:end end],it-1);
    X(:,it) = R(c+1);
    % display
    clf; 
    imagesc(1-X'); axis tight;  axis equal;  axis off;
    drawnow; 
    colormap gray(256);
    % save
    imwrite(Upsc(rescale(1-X'),2), [rep 'anim-' znum2str(it,3) '.png']);
end



