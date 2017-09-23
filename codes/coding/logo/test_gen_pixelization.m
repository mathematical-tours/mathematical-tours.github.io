%%
% Generate pixelization.

name = 'shannon';
rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end

addpath('./images/');
base_path = '/Users/gabrielpeyre/Dropbox/github/numerical-tours/matlab/';
addpath([base_path 'toolbox_general/']);
addpath([base_path 'toolbox_signal/']);

% target size
N = 512*2;

f = load_image(name);
r = size(f,1)/size(f,2);
P = round(N/r);
f = rescale(mean(f,3));
f = image_resize(f, [N P]);
f = rescale(mean(f,3));

Qlist = [40 20 10];
F = [f];
for i=1:length(Qlist)
    g = approx_constant(f,Qlist(i));
    F = [F,g];
end

imwrite(rescale(F), [rep name '-' num2str(N) '.png'], 'png');