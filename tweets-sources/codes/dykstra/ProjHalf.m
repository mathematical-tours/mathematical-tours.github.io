function x = ProjHalf(x,u)

if real(x*conj(u))>0
    x = x - real(x*conj(u))*u/abs(u)^2;
end

end