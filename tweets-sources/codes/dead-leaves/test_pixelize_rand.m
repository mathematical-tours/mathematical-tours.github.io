%%
% Generate pixelization.


name = 'paris2';
name = 'eye';
name = 'paris-cover'; 

addpath('../toolbox/');
rep = MkResRep(name);

addpath('./images/');
addpath('./toolbox/');

% target height
height = 2100;
height = 512/2;

remap = @(t)1.5 * ( (t(1)-1/2)^2 + (t(2)-1/2)^2 );
switch name
    case {'paris1' 'paris3'}
        remap = @(t)1.5 * (t(1)-.3)^2 + 3*(t(2)-1/2)^2;
    case {'paris6'}
        remap = @(t)1.5 * (t(1)-.3)^2 + 2*(t(2)-.1)^2;
    case 'paris4'
        remap = @(t)1.5 * (t(1)-.3)^2 + 3*(t(2)-.6)^2;
    case 'paris2'
        remap = @(t)1.5 * (t(1)-.5)^2 + 3*(t(2)-.35)^2;
    case 'gp'
        remap = @(t)min( 1*sqrt( (t(1)-.5)^2 + (t(2)-.62)^2 ), .4 ); 
    case 'eye'
        remap = @(t)3 * (t(1)-.36)^2 + 5*(t(2)-.54)^2;
    case 'paris-cover'
        remap = @(t)1 * abs(t(1)-.5)^2 + 3*abs(t(2)-.68)^2;
end

f = double(imread([name '.jpg']));
options.remap = remap;
options.nsamples = 50*10000;
options.nsamples = 3500*6;
options.remap = @(t).3;
options.remap = @(t).06+1 * abs(t(1)-.5)^2 + 3*abs(t(2)-.68)^2;
options.rep = rep;
[g,f0] = pixelize(f, height, options);

imageplot(g);
%imwrite(rescale(g), [rep name '-' num2str(height) '.png'], 'png');
% imwrite(rescale(f0), [rep name '-original.png'], 'png');
