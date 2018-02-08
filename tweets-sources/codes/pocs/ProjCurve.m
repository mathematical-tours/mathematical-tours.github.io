function y = ProjCurve(x,c)

[~,i] = min(abs(x-c));
y = c(i);

end