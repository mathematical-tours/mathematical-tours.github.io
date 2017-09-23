rep = 'results/quantif/';
if not(exist(rep))
    mkdir(rep);
end

addpath('toolbox/');

%%
% Parameters.

n = 256;
name = 'hibiscus';

%%
% Helpers.

Quant = @(x,q)min(floor( rescale(x)*q  ), q-1);


%%
% Plot the entropy function. 

t = linspace(1e-10,1,1000);
plot(t, -t.*log2(t))
axis tight; 
saveas(gcf, [rep 'entropy.eps'], 'epsc');

%%
% Load image.

f = load_image(name, n);
f = rescale(sum(f,3));
clf;
imageplot(f);

%%
% Zoom.

selx = 19:24; sely = 62:67;
selx = 154:154+4; sely = 118:118+4;
fz = f(selx,sely);

p = length(selx);
k = 32;
s = floor(1:1/k:p+1-1/k);
Upsc = @(x)x(s,s);

clf;
imageplot({f fz});

imwrite(f, [rep 'original.png'], 'png');
imwrite(Upsc(fz), [rep 'original-zoom.png'], 'png');

%%
% Quantify.

fprintf('--- original ---\n', q);
disp(fz);

for q=[2 4 8 16 255]
    a = Quant(f,q);
    az = Quant(fz,q);
    imwrite(a/(q-1), [rep 'q' num2str(q) '.png'], 'png');
    imwrite(Upsc(az)/(q-1), [rep 'q' num2str(q) '-zoom.png'], 'png');
    fprintf('--- q=%d ---\n', q);
    disp(az);
end
