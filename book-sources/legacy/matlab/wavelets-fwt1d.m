Jmax = log2(n)-1; Jmin = 0; fw = f;
for j=Jmax:-1:Jmin
    Coarse = subsampling(cconv(fw(1:2^(j+1)),h));
    Detail = subsampling(cconv(fw(1:2^(j+1)),g));
    fw(1:2^(j+1)) = cat(1, Coarse, Detail );
end