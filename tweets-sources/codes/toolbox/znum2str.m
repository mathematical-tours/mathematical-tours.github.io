function s = znum2str(x,k)

% force to happend trailing zeros to have string of length k

if x>=10^k
    warning('Problem with znum2str');
end

s = num2str(x);
m = k-1 - floor(log10(x));
for i=1:m
    s = ['0' s];
end
