M1 = MW;
for j=Jmin:Jmax
    A = M1(1:2^(j+1),1:2^(j+1));
    for d=1:2
        if d==1
            Coarse = A(1:2^j,:);
            Detail = A(2^j+1:2^(j+1),:);
        else
            Coarse = A(:,1:2^j);
            Detail = A(:,2^j+1:2^(j+1));                
        end
        Coarse = cconv(upsampling(Coarse,d),reverse(h),d);
        Detail = cconv(upsampling(Detail,d),reverse(g),d);
        A = Coarse + Detail;
    end
    M1(1:2^(j+1),1:2^(j+1)) = A;
end