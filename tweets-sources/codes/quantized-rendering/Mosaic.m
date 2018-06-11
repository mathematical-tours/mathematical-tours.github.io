%%
% You have to download the CIFAR-100 dataset from here
%   https://www.cs.toronto.edu/~kriz/cifar.html

addpath('../toolbox/');
rep = MkResRep();


if not(exist('Z'))
    Z = LoadCIFAR();
end

if 0
    % remove images with too large variance
    a = squeeze( std(std(std(Z,1,1),1,2),1,3) );
    [tmp,I] = sort(a);
    Z = Z(:,:,:,I(1:round(3*m/4)));
    m = size(Z,4);
end


% input image
repin = 'images/';
if not(exist('name'))
    name = 'cm';
end
f = load_image([repin name]);

% name = 'zoom';
% f = imread([name '.gif']);

% downsample input image
kx = 60; % number of horizontal squares
f1 = GenerateMosaic(f,Z,kx);

imageplot(f1);

% save to file
% repout = ['/Users/' username '/Dropbox/photos/2016/2016-04-08-mosaiques/'];
% imwrite(rescale(f1), [repout name '-mosaique.jpg'], 'jpg');
