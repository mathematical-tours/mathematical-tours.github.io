rep = '../results/quantization/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

rescale = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));

%%
% Parameters.

n = 256;
name = 'shannon';
f = load_image(name, n);
f = rescale(f);
imwrite(f, [rep 'original.png'], 'png');

%%
% Piecewise constant approximation.

for m=[8 16 32 64 128]
    imwrite(Approx(f,m), [rep 'approx-' num2str(m) '.png'], 'png');
end

%%
% Quantify.

Quant = @(x,q)min(floor( rescale(x)*q  ), q-1);
for q=[2 3 4 5 6 7 8 12 16 255]
    a = Quant(f,q);
    imwrite(a/(q-1), [rep 'quantize-' num2str(q) '.png'], 'png');
end
