Jmax = log2(n)-1; Jmin = 0;
MW = M;
for j=Jmax:-1:Jmin
    A = MW(1:2^(j+1),1:2^(j+1));
    for d=1:2
        Coarse = subsampling(cconv(A,h,d),d);
        Detail = subsampling(cconv(A,g,d),d);
        A = cat3(d, Coarse, Detail );
    end
    MW(1:2^(j+1),1:2^(j+1)) = A;
end