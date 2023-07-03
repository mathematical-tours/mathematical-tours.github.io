function z = FermatPoint(x)

z = mean(x);
for it=1:50
    w = 1./(1e-5 + abs(x-z)); w = w/sum(w);
    z = sum(w.*x);
end

end