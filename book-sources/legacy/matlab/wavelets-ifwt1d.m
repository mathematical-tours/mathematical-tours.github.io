f1 = fw;
for j=Jmin:Jmax
    Coarse = f1(1:2^j);
    Detail = f1(2^j+1:2^(j+1));
    Coarse = cconv(upsampling(Coarse,1),reverse(h),1);
    Detail = cconv(upsampling(Detail,1),reverse(g),1);
    f1(1:2^(j+1)) = Coarse + Detail;
end