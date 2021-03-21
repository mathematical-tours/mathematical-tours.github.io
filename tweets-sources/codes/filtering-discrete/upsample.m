function h1 = upsample(h,k)

h1 = zeros(size(h)*k);
h1(k:k:end,k:k:end) = h;

end