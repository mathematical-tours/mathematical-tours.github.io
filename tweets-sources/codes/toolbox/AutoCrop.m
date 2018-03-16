function AutoCrop(rep,name)

% Special cropping of images

basename = [rep name]; % 'interp-'
U = dir([basename '*.png']);

% load first image
f = imread([rep U(1).name]);
clf; hold on;
imagesc(f); axis image; axis off;
[a,b] = ginput(2);
a = sort(round(a)); b = sort(round(b));
for k=1:length(U)
    f = imread([rep U(k).name]);
    imwrite( f(b(1):b(2), a(1):a(2), :), [rep U(k).name], 'png');
end

end