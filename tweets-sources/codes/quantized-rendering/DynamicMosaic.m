%%
% Dynamical mosaic

name = 'blossom';

addpath('../toolbox/');
rep = MkResRep(name);

if not(exist('Z'))
    Z = LoadCIFAR();
    Z = Z/255;
end

[F,cm] = imread(['images/' name '.gif']);

q = min(size(F,4),80);
kx = 30; % number of horizontal squares
    
for i=1:q
    progressbar(i,q);
    f = cm(F(:,:,:,i)+1,:);
    f = reshape(f, [size(F,1) size(F,2) 3]);
    % downsample input image
    f1 = GenerateMosaic(f,Z,kx);
    clf; imageplot(f1); drawnow;
    % save
    imwrite(f, [rep 'original-' znum2str(i,3) '.png'], 'png');
    imwrite(f1, [rep 'mosaic-' znum2str(i,3) '.png'], 'png');
end