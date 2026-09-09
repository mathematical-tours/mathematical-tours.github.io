function y = decode_rm(x)
N = length(x);
f = fwt((-1).^x);
[v,a] = max(abs(f));
y = encode_rm( a-1,f(a)<0,log2(N) );