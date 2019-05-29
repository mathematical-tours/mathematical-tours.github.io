name = 'hibiscus';

addpath('../toolbox/');
rep = MkResRep();

n = 256;
f = load_image(name,n);
f = clamp(f/255);

imwrite( f, [rep 'original.png'] );

imwrite( clamp(f+.1*randn(n,n,3)), [rep 'noisy.png'] );