%%
% Color quantization.

name = 'landscape';
addpath('../toolbox/');
rep = MkResRep(name);

n = 512; 

f = load_image(name, n);
f = rescale(f);
imwrite(f, [rep 'original.png'], 'png');

q = 50
mlist = round(1+200*linspace(0,1,q).^1.5);
for m=1:q
    g = Approx(f,mlist(m));
    clf; imageplot(g); drawnow;
    imwrite(g, [rep 'approx-' znum2str(m,3) '.png'], 'png');
end


%%
% Quantify.

for m=1:q
    g = QuantizeColor(f,m);
    clf; imageplot(g); drawnow;
    imwrite(g, [rep 'quantize-' znum2str(m,3) '.png'], 'png');
end