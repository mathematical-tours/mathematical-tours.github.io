function y = encode_rm(a,b,n)
y = zeros(2^n,1); y(a+1) = 1;
y = fwt(y); y = (1-y)/2;
if(b==1) y = 1-y; end;