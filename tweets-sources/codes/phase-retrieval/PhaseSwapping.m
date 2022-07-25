addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);

n = 256*2;


% windowing
t = linspace(0,1,n)';
t = sin(pi*t).^2;

name = {'shannon' 'hibiscus'};
f = {}; F = {};
for k=1:2
    f{k} = rescale( sum(load_image(name{k}, n),3) );
    f{k} = f{k}.*(t*t');
    F{k} = fft2(f{k});
end

PhSwap = @(a,b) real( ifft2( abs(a).*b./abs(b)  ) );
g = { PhSwap(F{1},F{2}), PhSwap(F{2},F{1})  };

p = round(n*.8);
f = crop(f,p);
g = crop(g,p);

clf;
imageplot({f{:},g{:}},{'f', 'g', 'abs(f) phase(g)', 'abs(g) phase(f)'}, 2,2);

for k=1:2
    imwrite(clamp(f{k}), [rep 'input-' num2str(k) '.png']);
    imwrite(clamp(g{k}), [rep 'swap-' num2str(k) '.png']);
end