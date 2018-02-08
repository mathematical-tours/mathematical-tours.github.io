rep = '../results/hist-eq/';
[~,~] = mkdir(rep);

rescale = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)));
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

name = 'hibiscus.png';
f = imread(name);
f = rescale( mean(f,3) ).^.5;
f0 = f;



[fS,I] = sort(f0(:)); f(I) = linspace(0,1,length(f(:)));

imagesc(f0);
colormap gray(256);
imagesc(f);
colormap gray(256);

imwrite(f0, [rep 'original.png']);

q = 9;
for k=1:q
    t = (k-1)/(q-1);
    f(I) = (1-t)*fS + t*linspace(0,1,length(f(:)))';
    imwrite(f, [rep 'equalized-' num2str(k) '.png']);
    % disp
    imagesc(f);
    colormap gray(256);
    drawnow;
    % histogram
    m = 25;
    h = hist(f, m); h = h/sum(h);
    if k==1
        vmax=max(h);
    end
    clf;
    bar(h, 'FaceColor', [t 0 1-t]);
    axis([1 m 0 vmax*1.01]);
    SetAR(1/3); 
    set(gca, 'XTick', [], 'YTick', []);
    saveas(gcf, [rep 'equalize-' num2str(k) '.eps'], 'epsc');
end