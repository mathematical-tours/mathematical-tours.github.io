% test des DCT

n = 16;
x = rand(n,1);

y1 = dct2(x);
y2 = dct2_triv(x);
y3 = dct3(y1);
y4 = dct3_triv(x);
