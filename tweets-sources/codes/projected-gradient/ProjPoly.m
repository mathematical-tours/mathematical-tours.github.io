function z = ProjPoly(z,Cc,E)

if inpolygon(imag(z),real(z),real(E),imag(E))
    % nothing
else
    [~,i] = min(abs(z-Cc));
    z = Cc(i);
end

end