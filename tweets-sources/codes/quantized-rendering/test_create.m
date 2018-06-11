username = 'gabrielpeyre'; 
username = 'gpeyre';

% load library;
if not(exist('Z'))
    proot = ['/Users/' username '/Documents/MATLAB/data/'];    
    load([proot 'cifar-100-matlab/train']);
    a = data;
    load([proot 'cifar-100-matlab/test']);
    a = [a; data];
    Z = reshape(a', [32 32 3 size(a,1)]);
    Z = double(Z);
    clear data; clear a;
end

n = size(Z,1);
m = size(Z,4);

addpath('toolbox/');

if 0
    % remove images with too large variance
    a = squeeze( std(std(std(Z,1,1),1,2),1,3) );
    [tmp,I] = sort(a);
    Z = Z(:,:,:,I(1:round(3*m/4)));
    m = size(Z,4);
end

% means
mu = mean(mean(Z,1),2);
mu = squeeze(mu)';

repin = 'images/';
if not(exist('name'))    
    name = '2015-12-25-brest-023';
end
f = load_image([repin name]);

kx = 250; % number of horizontal squares
w = floor(size(f,1)/kx);
ky = floor(size(f,2)/w);
f = f(1:w*kx,1:w*ky,:,:);
k = kx*ky;

g = reshape(f, [w kx w ky 3]);
g = squeeze( mean(mean(g,1),3) );

% find the best index
G = reshape(g, [k 3]);

% randomization factor
q = 400; 

if 0
    D = distmat(G',mu');
    [tmp,I] = sort(D,2);
else
    I = knnsearch(G,mu,q);
end


% randomize a little bit to avoid repetition

ind = floor(rand(k,1)*q)+1;
J = [];
for i=1:k
    J(i) = I(i,ind(i));
end
J = reshape(J,[kx ky]);

% re-create the final image
f1 = zeros(kx*n, ky*n,3);
for i=1:kx
    progressbar(i,kx);
    for j=1:ky
        q = J(i,j); % selected patch
        x = (i-1)*n+1:i*n;
        y = (j-1)*n+1:j*n;
        a = double( Z(:,:,:,q) );
        a = permute(a, [2 1 3]);
        a = a + repmat( reshape(-mu(q,:),[1 1 3]) + g(i,j,:) , [n n] );
        a = clamp(a,0,255);
        f1(x,y,:) = a;
    end
end
repout = ['/Users/' username '/Dropbox/photos/2016/2016-04-08-mosaiques/'];
imwrite(rescale(f1), [repout name '-mosaique.jpg'], 'jpg');
