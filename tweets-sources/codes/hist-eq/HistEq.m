
addpath('../toolbox/');
rep = MkResRep();


rescale = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

name = 'hibiscus.png';
f = imread(name);
f = rescale( mean(f,3) ).^1.5;
f0 = f;


[fS,I] = sort(f0(:)); f(I) = linspace(0,1,length(f(:)));

imagesc(f0);
colormap gray(256);
imagesc(f);
colormap gray(256);

imwrite(f0, [rep 'original.png']);

% target histogram
G = clamp( sort( .8 + .2*randn(length(f(:)),1)*.2 ), 0,1);
G = linspace(0,1,length(f(:)))';


q = 50;
for k=1:q
    t = (k-1)/(q-1);
    f(I) = (1-t)*fS + t*G;
    imwrite(f, [rep 'img-' znum2str(k,2) '.png']);
    % disp
    imagesc(f);
    colormap gray(256);
    drawnow;
    % histogram
    m = 40;
    x = linspace(0,1,m);
    h = hist(f(:), x); h = h/sum(h);
    if k==1
        vmax=max(h);
    end
    
    clf;
    bar(x, h, 'FaceColor', [t 0 1-t]);
    axis([0 1 0 vmax*1.01]);
    SetAR(1/3); 
    set(gca, 'XTick', [], 'YTick', []);
    saveas(gcf, [rep 'hist-' znum2str(k,2) '.png'], 'png');
end

% AutoCrop(rep, ['hist-']); 